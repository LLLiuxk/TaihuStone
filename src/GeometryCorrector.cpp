// TopoOpt.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#include <string>
#include <fstream>
#include <iomanip>
#include "GeometryCorrector.h"



GeometryCorrector::GeometryCorrector(int logical_w, int logical_h, int pad, double filter_radius)
    : nx_logical(logical_w), ny_logical(logical_h), padding(pad), r_min(filter_radius) {

    nx = nx_logical + 2 * padding;
    ny = ny_logical + 2 * padding;

    int size = nx * ny;
    rho.resize(size, 0.0);
    rho_init.resize(size, 0.0);
    sensitivity.resize(size, 0.0);
    x.resize(size, 0.0);         // 初始设为 0.5
    rho_tilde.resize(size, 0.0);

    beta = 0.5;                  // 初始 beta 较小，保持线性
    eta = 0.5;
}

    // 初始化：创建一个简单的测试形状（例如一个悬空的圆，不可打印）
void GeometryCorrector::initializeTestShape() {
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            // 在中心画一个圆
            double dx = x - nx / 2.0;
            double dy = y - ny / 2.0; // 假设 Y 轴向上
            if (dx * dx + dy * dy < (nx / 4.0) * (nx / 4.0)) {
                setDensity(rho_init, x, y, 1.0);
                setDensity(rho, x, y, 1.0); // 初始猜测 = 目标形状
            }
            else {
                setDensity(rho_init, x, y, 0.0);
                setDensity(rho, x, y, 0.0);
            }
        }
    }
    std::cout << "Initialized test shape (Circle)." << std::endl;
    x = rho_init;
}

void GeometryCorrector::initTShape() {

    std::fill(x.begin(), x.end(), 0.0);
    std::fill(rho.begin(), rho.end(), 0.0);
    std::fill(rho_init.begin(), rho_init.end(), 0.0);

    for (int y = ny - padding; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            setDensity(this->x, x, y, 0.0);        // 设计变量设为1
            setDensity(rho, x, y, 0.0);      // 物理密度设为1
            setDensity(rho_init, x, y, 0.0); // 目标也设为1 (避免目标函数报错)
        }
    }

    int mid_x = nx_logical / 2;
    int beam_width = nx_logical / 5;
    int top_height = ny_logical / 5;

    for (int j = 0; j < ny_logical; ++j) {
        for (int i = 0; i < nx_logical; ++i) {
            // 逻辑坐标转物理坐标
            int phys_x = i + padding;
            int phys_y = j + padding;

            bool is_vertical_post = (i >= mid_x - beam_width / 2 && i < mid_x + beam_width / 2) && (j >= top_height);
            bool is_horizontal_bar = (j < top_height);

            if (is_vertical_post || is_horizontal_bar) {
                setDensity(x, phys_x, phys_y, 1.0);
                setDensity(rho, phys_x, phys_y, 1.0);
                setDensity(rho_init, phys_x, phys_y, 1.0);
            }
        }
    }

    std::cout << "Initialized test shape (T)." << std::endl;

}

// 辅助函数：设置密度（带边界检查）
void GeometryCorrector::setDensity(std::vector<double>& field, int x, int y, double val) {
    if (x >= 0 && x < nx && y >= 0 && y < ny) {
        field[y * nx + x] = val;
    }
}

// 辅助函数：获取密度
double GeometryCorrector::getDensity(const std::vector<double>& field, int x, int y) {
    if (x >= 0 && x < nx && y >= 0 && y < ny) {
        return field[y * nx + x];
    }
    return 0.0; // 边界外视为无材料
}

// ---------------------------------------------------------
// 核心逻辑 1: 计算目标函数 (几何偏差) 及其导数
// J = Sum( (rho_i - rho_init_i)^2 )
// dJ/drho = 2 * (rho_i - rho_init_i)
// ---------------------------------------------------------
double GeometryCorrector::computeObjectiveAndGradient() {
    //double total_diff = 0.0;

    //for (int i = 0; i < nx * ny; ++i) {
    //    double diff = rho[i] - rho_init[i];

    //    // 目标函数值累加
    //    total_diff += diff * diff;

    //    // 计算梯度 (这是最内层的导数，之后需要乘上滤波和投影的链式法则)
    //    sensitivity[i] = 2.0 * diff;
    //}

    //return total_diff;

    double total_cost = 0.0;

    // --------------------------------------------------------
    // 参数调整区
    // --------------------------------------------------------

    // 1. 解决“削减过度”：调大保真权重
    // 之前可能太小了，导致T型臂稍有悬垂就被删光。
    // 建议：20.0 -> 50.0 左右尝试
    double match_weight = Match_weight;

    // 2. 基础体积费：原本的微小生存压力
    double vol_weight = Vol_weight;

    // 3. [新增] 体积上限惩罚
    // 允许的最大体积占比 (例如 40%)
    double max_vol_frac = Vol_frac;
    // 惩罚力度 (要非常大，像一堵墙)
    double vol_barrier_weight = 50.0;

    // --------------------------------------------------------

    // 先计算当前体积
    double current_vol_sum = 0.0;
    for (double v : rho) current_vol_sum += v;
    double current_vol_frac = current_vol_sum / (nx * ny);

    // 计算体积超标带来的额外梯度 (常数项)
    double vol_barrier_grad = 0.0;
    if (current_vol_frac > max_vol_frac) {
        // 这是一个二次惩罚函数： Cost += W * (V - V_max)^2
        double diff = current_vol_frac - max_vol_frac;
        total_cost += vol_barrier_weight * diff * diff * (nx * ny);

        // 导数： dC/dx = 2 * W * (V - V_max) * (1/N)
        // 这个梯度是正的（正值代表由于体积太大，每个单元都要受罚）
        vol_barrier_grad = 2.0 * vol_barrier_weight * diff;// / (nx * ny);
    }

    // 逐单元计算梯度
    for (int i = 0; i < nx * ny; ++i) {
        double r = rho[i];
        double r0 = rho_init[i];

        // --- A. 保真项 (One-sided) ---
        // 只有当 r < r0 (丢失材料) 时才惩罚
        double diff_match = std::max(0.0, r0 - r);
        total_cost += match_weight * diff_match * diff_match;
        double grad_match = (r < r0) ? (-2.0 * match_weight * diff_match) : 0.0;

        // --- B. 基础体积项 ---
        total_cost += vol_weight * r;
        double grad_base_vol = vol_weight;

        // --- 总梯度 ---
        // 梯度 = 保真 + 基础体积 + 体积上限惩罚
        sensitivity[i] = grad_match + grad_base_vol + vol_barrier_grad;
    }

    // 用于Debug：看看到底是不是体积超标了
    // std::cout << "Vol: " << current_vol_frac << " (Max: " << max_vol_frac << ")" << std::endl;

    return total_cost;
}

// ---------------------------------------------------------
// 核心逻辑 2: 密度滤波 (Density Filter)
// 作用：平滑设计变量，引入长度尺度，避免棋盘格。
// 公式：rho_tilde_i = Sum(w_ij * x_j) / Sum(w_ij)
// ---------------------------------------------------------
void GeometryCorrector::applyFilter() {
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            double sum_val = 0.0;
            double sum_weight = 0.0;
            int idx_current = i * nx + j;
            int r_ceil = static_cast<int>(std::ceil(r_min));

            for (int dy = -r_ceil; dy <= r_ceil; ++dy) {
                for (int dx = -r_ceil; dx <= r_ceil; ++dx) {
                    int ni = i + dy;
                    int nj = j + dx;
                    // 这里的边界检查是针对“物理边界”的
                    // 由于有 padding，逻辑部分的边界现在处于物理中心安全区
                    if (ni >= 0 && ni < ny && nj >= 0 && nj < nx) {
                        double dist = std::sqrt(dx * dx + dy * dy);
                        if (dist < r_min) {
                            double weight = r_min - dist;
                            sum_val += weight * x[ni * nx + nj];
                            sum_weight += weight;
                        }
                    }
                }
            }
            rho_tilde[idx_current] = sum_val / (sum_weight + 1e-10);
        }
    }
}

// ---------------------------------------------------------
// 核心逻辑 3: Heaviside 投影 (Projection)
// 作用：将灰度推向 0 或 1。
// 公式：使用 tanh 近似 Heaviside 函数
// ---------------------------------------------------------
void GeometryCorrector::applyProjection() {
    for (int i = 0; i < nx * ny; ++i) {
        double num = std::tanh(beta * eta) + std::tanh(beta * (rho_tilde[i] - eta));
        double den = std::tanh(beta * eta) + std::tanh(beta * (1.0 - eta));
        rho[i] = num / den;
    }
}

// ---------------------------------------------------------
// 核心逻辑 4: 链式法则求导 (Chain Rule)
// 目标：计算 dJ/dx (对设计变量的导数)，用于梯度下降。 
// 路径：dJ/dx = (dJ/drho) * (drho/drho_tilde) * (drho_tilde/dx)
// ---------------------------------------------------------
void GeometryCorrector::applyChainRule() {
    std::vector<double> sensitivity_tilde(nx * ny, 0.0);

    // 1. dJ/drho -> dJ/drho_tilde
    for (int i = 0; i < nx * ny; ++i) {
        double num_grad = beta * (1.0 - std::pow(std::tanh(beta * (rho_tilde[i] - eta)), 2));
        double den = std::tanh(beta * eta) + std::tanh(beta * (1.0 - eta));
        double d_projection = num_grad / den;
        sensitivity_tilde[i] = sensitivity[i] * d_projection;
    }

    // 2. dJ/drho_tilde -> dJ/dx (反向滤波)
    std::fill(sensitivity.begin(), sensitivity.end(), 0.0);
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            int idx_tilde = i * nx + j;

            // 重新计算权重和 (归一化因子)
            double sum_weight = 0.0;
            int r_ceil = static_cast<int>(std::ceil(r_min));
            for (int dy = -r_ceil; dy <= r_ceil; ++dy) {
                for (int dx = -r_ceil; dx <= r_ceil; ++dx) {
                    int ni = i + dy, nj = j + dx;
                    if (ni >= 0 && ni < ny && nj >= 0 && nj < nx) {
                        double dist = std::sqrt(dx * dx + dy * dy);
                        if (dist < r_min) sum_weight += (r_min - dist);
                    }
                }
            }

            double val_to_distribute = sensitivity_tilde[idx_tilde] / (sum_weight + 1e-10);

            // 分发梯度
            for (int dy = -r_ceil; dy <= r_ceil; ++dy) {
                for (int dx = -r_ceil; dx <= r_ceil; ++dx) {
                    int ni = i + dy, nj = j + dx;
                    if (ni >= 0 && ni < ny && nj >= 0 && nj < nx) {
                        double dist = std::sqrt(dx * dx + dy * dy);
                        if (dist < r_min) {
                            double weight = r_min - dist;
                            int idx_x = ni * nx + nj;
                            sensitivity[idx_x] += val_to_distribute * weight;
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------
// 核心逻辑 5: 悬垂约束 (Overhang Constraint)
// 这是一个简化的几何约束：惩罚那些“上方有材料但下方无支撑”的单元。
// ---------------------------------------------------------
//void GeometryCorrector::applyOverhangConstraint(double penalty_weight) {
//    // 遍历物理区域 (除最后一行)
//    for (int i = 0; i < ny - 1; ++i) {
//        for (int j = 0; j < nx; ++j) {
//            int idx_curr = i * nx + j;
//            int idx_below = (i + 1) * nx + j;
//
//            double r_c = rho[idx_curr];
//            double r_b = rho[idx_below];
//
//            // 这里的逻辑现在很安全：
//            // 如果结构长到了倒数第二行 (i=ny-2)，它的下方是最后一行 (i=ny-1)。
//            // 因为我们在 initTShapeWithPadding 把最后一行设为了 1.0 (基板)，
//            // 所以 (1 - r_b) 为 0，悬垂惩罚消失。
//            // 这允许结构安全地“着陆”。
//
//            double dP_drc = penalty_weight * std::pow(1.0 - r_b, 2);
//            double dP_drb = -2.0 * penalty_weight * r_c * (1.0 - r_b);
//
//            sensitivity[idx_curr] += dP_drc;
//            sensitivity[idx_below] += dP_drb;
//        }
//    }
//}

void GeometryCorrector::applyOverhangConstraint(double penalty_weight) {
    // 遍历除最后一行外的所有点
    for (int i = 0; i < ny - 1; ++i) {
        for (int j = 0; j < nx; ++j) {
            int idx_curr = i * nx + j;
            double r_c = rho[idx_curr];

            // --- 关键修改：寻找最强的支撑者 ---

            // 1. 获取下方三个邻居的密度
            // 左下方 (Bottom-Left)
            double r_bl = (j > 0) ? rho[(i + 1) * nx + (j - 1)] : 0.0;
            // 正下方 (Bottom-Center)
            double r_bc = rho[(i + 1) * nx + j];
            // 右下方 (Bottom-Right)
            double r_br = (j < nx - 1) ? rho[(i + 1) * nx + (j + 1)] : 0.0;

            // 2. 计算支撑度 Support Value (S)
            // S = max(r_bl, r_bc, r_br)
            // 意思就是：只要这三个里有一个是硬的，我就能踩着它往上长
            double support_val = r_bc;
            int support_idx = (i + 1) * nx + j; // 记录谁是“大腿”

            if (r_bl > support_val) {
                support_val = r_bl;
                support_idx = (i + 1) * nx + (j - 1);
            }
            if (r_br > support_val) {
                support_val = r_br;
                support_idx = (i + 1) * nx + (j + 1);
            }

            // 3. 计算惩罚 P = w * r_c * (1 - S)^2
            // 物理意义：如果我有材料(r_c=1)，且最强的支撑者也很弱(S=0)，则惩罚巨大
            double term = 1.0 - support_val;
            double P = penalty_weight * r_c * term * term; // 这里只是公式演示，实际不需要加到 sensitivity

            // 4. 计算梯度 (Sensitivity)

            // A. 对当前点 r_c 的导数: dP/dr_c = w * (1 - S)^2
            // 正值：如果下方支撑不够，我就该被删掉
            double dP_drc = penalty_weight * term * term;
            sensitivity[idx_curr] += dP_drc;

            // B. 对支撑点 r_support 的导数: dP/dS = -2 * w * r_c * (1 - S)
            // 负值：如果我需要支撑，支撑点应该增加材料
            // 注意：梯度只传给那个“最大值”的邻居（Max 函数的导数特性）
            // 这会引导优化器去加强那个已经是最大的邻居，形成明确的生长路径
            double dP_dS = -2.0 * penalty_weight * r_c * term;

            // 这里的判断是为了处理边界 padding 或基板 (假设基板不参与优化)
            if (support_idx < nx * ny) {
                sensitivity[support_idx] += dP_dS;
            }
        }
    }
}

void GeometryCorrector::resetBoundarySensitivities() {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            // 判断是否在 Padding 区域
            // 逻辑区域是 [padding, nx-padding)
            bool is_padding = (i < padding) || (i >= nx - padding) ||
                (j < padding) || (j >= ny - padding);

            if (is_padding) {
                // 强制梯度为 0，锁定该变量
                sensitivity[j * nx + i] = 0.0;

                // 可选：同时也强制把 x 的值重置回去，防止数值漂移
                // 但通常梯度为0就足够了
            }
        }
    }
}

void GeometryCorrector::printTopology() {
    std::cout << "\n--- Current Topology (Beta=" << beta << ") ---\n";
    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            double val = rho[i * nx + j];
            if (val > 0.75) std::cout << "#";      // 实体
            else if (val > 0.25) std::cout << "+"; // 灰色过渡带
            else std::cout << ".";                 // 空洞
        }
        std::cout << "\n";
    }
    std::cout << "-----------------------------------\n";
}

// ---------------------------------------------------------
// 工具：导出为 PPM 图片 (无需库即可查看)
// ---------------------------------------------------------
void GeometryCorrector::exportPPM(const std::string& filename, const std::vector<double>& field) {
    std::ofstream ofs(filename, std::ios::binary);
    ofs << "P6\n" << nx << " " << ny << "\n255\n";
    for (int i = 0; i < nx * ny; ++i) {
        unsigned char val = static_cast<unsigned char>((1.0 - field[i]) * 255); // 0=黑(材料), 1=白(空)
        ofs << val << val << val;
    }
    ofs.close();
    // std::cout << "Saved " << filename << std::endl;
}

void GeometryCorrector::imposeSymmetry() {
    // 假设设计域关于 x 轴中心对称
    // 注意：这里我们只处理物理区域 [0, nx)
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx / 2; ++x) { // 只遍历左半边
            int idx_left = y * nx + x;
            int idx_right = y * nx + (nx - 1 - x); // 对应的右半边位置

            // 取平均值
            double avg_val = (this->x[idx_left] + this->x[idx_right]) / 2.0;

            // 强行赋值给两边
            this->x[idx_left] = avg_val;
            this->x[idx_right] = avg_val;
        }
    }
}
   
void gradientDescentUpdate(std::vector<double>& x,
    const std::vector<double>& dc,
    double learning_rate) {
    //for (size_t i = 0; i < x.size(); ++i) {
    //    // x_new = x - learning_rate * sensitivity
    //    // 因为我们定义的是“代价函数”的梯度，所以我们要沿着负梯度方向走
    //    x[i] = x[i] - learning_rate * dc[i];

    //    // 投影到 [0, 1]
    //    if (x[i] < 0.0) x[i] = 0.0;
    //    if (x[i] > 1.0) x[i] = 1.0;
    //}
    double damping = 0.5;

    for (size_t i = 0; i < x.size(); ++i) {
        // 计算原本计划的更新值
        double proposed_val = x[i] - learning_rate * dc[i];

        // 截断
        if (proposed_val < 0.0) proposed_val = 0.0;
        if (proposed_val > 1.0) proposed_val = 1.0;

        // --- 应用阻尼 (平滑更新) ---
        // 这一步能极大地抑制震荡
        x[i] = damping * x[i] + (1.0 - damping) * proposed_val;
    }
}

void optimalityCriteriaUpdate(std::vector<double>& x,
    const std::vector<double>& dc,
    double target_vol_frac) {
    //double l1 = 0.0, l2 = 100000.0, move = 0.2;
    //size_t n = x.size();
    //std::vector<double> xnew(n);

    //// 二分法寻找拉格朗日乘子 (Lagrange Multiplier)
    //while ((l2 - l1) / (l1 + l2 + 1e-10) > 1e-3) {
    //    double lmid = 0.5 * (l2 + l1);
    //    double current_vol = 0.0;

    //    for (size_t i = 0; i < n; ++i) {
    //        // OC 更新公式 (启发式): x_new = x * sqrt(-sensitivity / lambda)
    //        // 注意：因为我们是最小化目标，sensitivity 通常是负的。
    //        // 这里为了演示，假设 sensitivity 是正数代表“需要增加材料”
    //        // 严谨的物理公式中，柔度导数是负的，所以通常用 -dc[i]

    //        double term = x[i] * std::sqrt(std::max(1e-10, dc[i]) / lmid);

    //        // 移动限制 (Move Limit) 防止震荡
    //        double val = std::max(0.0, std::max(x[i] - move, std::min(1.0, std::min(x[i] + move, term))));
    //        xnew[i] = val;
    //        current_vol += val;
    //    }

    //    if (current_vol / n > target_vol_frac) {
    //        l1 = lmid; // 体积太大，增大惩罚
    //    }
    //    else {
    //        l2 = lmid; // 体积太小，减小惩罚
    //    }
    //}
    //x = xnew; // 更新设计变量
    // 
    // 学习率 (步长)，可能需要调试，先设小一点
    double learning_rate = 0.05;

    // 1. 梯度下降 (注意：我们要最小化 J，所以是 x = x - learning_rate * gradient)
    // 你的 dc[i] 是正值（代表代价增加），所以我们要减去它。
    std::vector<double> xnew = x;
    for (size_t i = 0; i < x.size(); ++i) {
        xnew[i] = x[i] - learning_rate * dc[i];

        // 简单截断 [0, 1]
        if (xnew[i] < 0.0) xnew[i] = 0.0;
        if (xnew[i] > 1.0) xnew[i] = 1.0;
    }

    // 2. 满足体积约束 (简单的二分法或者是单纯的归一化)
    // 这里做一个简单的整体平移/缩放来满足体积，或者更严谨地用二分法找由于投影引起的截距
    // 为了最简单的 Debug，我们先忽略严格的体积约束，只看形状变化：
    // (如果形状变了，说明梯度方向对了)

    // 稍微高级一点的体积修正（二分法找平移量 mu）：
    double l1 = -1.0, l2 = 1.0;
    while (l2 - l1 > 1e-4) {
        double mu = 0.5 * (l1 + l2);
        double current_vol = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            double val = x[i] - learning_rate * dc[i] + mu; // 引入平移量
            val = std::max(0.0, std::min(1.0, val));
            current_vol += val;
        }
        if (current_vol / x.size() > target_vol_frac) l2 = mu; // 体积太大，mu要小点
        else l1 = mu;
    }

    // 应用找到的 mu
    double mu = l1;
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = std::max(0.0, std::min(1.0, x[i] - learning_rate * dc[i] + mu));
    }
}