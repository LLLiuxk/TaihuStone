#include "TopoOpt.h"

void initTShape() {

    int nx_logical = 40;
    int ny_logical = 20;
    int padding = 0;
    int nx = nx_logical + 2 * padding;
    int ny = ny_logical + 2 * padding;

    std::vector<double> rho;
    int size = nx * ny;
    rho.resize(size, 0.0);
    //for (int y = ny - padding; y < ny; ++y) {
    //    for (int x = 0; x < nx; ++x) {
    //        if (x >= 0 && x < nx && y >= 0 && y < ny) {
    //            rho[y * nx + x] = 1.0;
    //        }
    //    }
    //}

    int mid_x = nx_logical / 2;
    int beam_width = nx_logical / 5;
    int top_height = ny_logical / 5;

    for (int j = 0; j < ny_logical + padding; ++j) {
        for (int i = 0; i < nx_logical; ++i) {
            // 逻辑坐标转物理坐标
            int phys_x = i + padding;
            int phys_y = j + padding;

            bool is_vertical_post = (i >= mid_x - beam_width / 2 && i < mid_x + beam_width / 2) && (j >= top_height);
            bool is_horizontal_bar = (j < top_height);

            if (is_vertical_post || is_horizontal_bar) {
                rho[phys_y * nx + phys_x] = 1.0;
            }
        }
    }

    std::cout << "Initialized test shape (T)." << std::endl;

    cv::Mat topology_img(ny, nx, CV_8UC1);

    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            double val = rho[i * nx + j];
            // 二值化：大于0.5为白色(255)，小于等于0.5为黑色(0)
            topology_img.at<uchar>(i, j) = (val > 0.5) ? 255 : 0;
        }
    }

    // 保存为JPG
    std::string filename = "topology.jpg";
    cv::imwrite(filename, topology_img);
    std::cout << "Topology image saved as: " << filename << std::endl;

}

// 在 TopologyOptimizer 类中添加
Eigen::VectorXd flatten(const Eigen::MatrixXd& mat) {
    // 将 MatrixXd (Rows x Cols) 展平为 VectorXd (Rows*Cols x 1)
    // Eigen 默认列优先，但在拓扑优化中，我们通常按行或列统一即可。
    // 这里我们将矩阵重塑为向量。
    Eigen::VectorXd vec(mat.size());
    // 按行优先展平 (为了和循环逻辑对应)
    int idx = 0;
    for (int i = 0; i < mat.rows(); ++i) {
        for (int j = 0; j < mat.cols(); ++j) {
            vec(idx++) = mat(i, j);
        }
    }
    return vec;
}

// 将 1D 向量还原为 2D 矩阵
void unflatten(const Eigen::VectorXd& vec, Eigen::MatrixXd& mat) {
    int idx = 0;
    for (int i = 0; i < mat.rows(); ++i) {
        for (int j = 0; j < mat.cols(); ++j) {
            mat(i, j) = vec(idx++);
        }
    }
}

// ==========================================
// 2. 拓扑优化核心类
// ==========================================


TopologyOptimizer::TopologyOptimizer(const std::string& image_path, Config config) : cfg(config) {
    // 读取图像并归一化
    cv::Mat img = cv::imread(image_path, cv::IMREAD_GRAYSCALE);
    if (img.empty()) throw std::runtime_error("Cannot load image");

    rows = img.rows;
    cols = img.cols;
    rho.resize(rows, cols);
    target_shape.resize(rows, cols);
    velocity = Eigen::MatrixXd::Zero(rows, cols); // 初始化为0

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = img.at<uint8_t>(i, j) / 255.0;
            rho(i, j) = val; // 初始设计变量
            target_shape(i, j) = val; // 目标形状
        }
    }
    rho_tilde = rho; // 初始占位

    // 初始化四元素模式矩阵 (Ref: Eq. 25 & 28)
    // B_l 用于左下邻域
    B_l << 0.5, -0.5, -0.5, 0.5,
        -0.5, -0.5, 0.5, 0.5,
        0.75, 0.25, -0.25, 0.25;

    // B_r 用于右下邻域
    B_r << -0.5, -0.5, 0.5, 0.5,
        -0.5, 0.5, 0.5, -0.5,
        0.5, 0.25, -0.25, 0.25;
}

void TopologyOptimizer::solve() {
    // 1. 初始化参数
    int n = rows * cols; // 变量数量
    int m = 3;           // 约束数量 (P 和 Q)
    beta = cfg.beta_start;

    // 2. 初始化 MMA 求解器
    // MMASolver(n_vars, n_constraints)
    MMASolver mma(n, m);

    // 3. 准备数据容器
    Eigen::VectorXd x(n);           // 设计变量 (当前)
    Eigen::VectorXd xold(n);        // 设计变量 (上一步)
    Eigen::VectorXd df0dx(n);       // 目标函数梯度
    Eigen::VectorXd constraints(m); // 约束值 [P, Q]

    // 约束梯度矩阵: m 行 n 列 (Row-Major flattening for MMA usually)
    // 或者是 n 行 m 列，取决于库的实现。
    // jdumas/mma 的 Update 接受 double* dgdx，通常是 m*n 大小的数组。
    // 它的排列顺序通常是: [df1/dx1, df1/dx2... | df2/dx1, df2/dx2...] (Constraint-Major)
    Eigen::MatrixXd dfdx(m, n);

    // 变量上下界
    Eigen::VectorXd xmin = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd xmax = Eigen::VectorXd::Ones(n);

    // 初始化 x
    x = flatten(rho);
    xold = x;

    std::cout << "Starting MMA Optimization..." << std::endl;
    printTopology();
    // 4. 优化主循环
    for (int iter = 0; iter < cfg.max_iter; ++iter) {

        // --- A. 更新物理场 ---
        unflatten(x, rho); // 将 MMA 的 x 同步回 rho

        // 施加对称性约束 (可选，但在每次更新变量后立即做最好)
        // imposeSymmetry(); // 如果你有这个函数
        
        Eigen::MatrixXd rho_bar = apply_linear_filter(rho, cfg.r_min);
        apply_projection(rho_bar, beta, cfg.eta); // 更新 rho_tilde

        // --- B. 计算目标函数及梯度 (J, dJ/dx) ---
        // 使用非对称目标函数逻辑
        Eigen::MatrixXd sens_J_mat = Eigen::MatrixXd::Zero(rows, cols);
        double J_val = 0.0;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double diff = rho_tilde(i, j) - target_shape(i, j);
                J_val += 0.5 * diff * diff; // 仅用于打印

                // 灵敏度计算
                if (diff < 0) sens_J_mat(i, j) = diff * cfg.cost_remove;
                else          sens_J_mat(i, j) = diff * cfg.cost_add;
            }
        }
        // 链式法则：dJ/drho
        Eigen::MatrixXd grad_J_mat = backprop_filter(sens_J_mat, rho_bar, beta);
        df0dx = flatten(grad_J_mat);

        // --- C. 计算约束及梯度 (P, Q, dP/dx, dQ/dx) ---

        // 约束 1: 悬垂角 P
        Eigen::MatrixXd sens_P_mat = Eigen::MatrixXd::Zero(rows, cols);
        double P_val = compute_overhang_constraint(sens_P_mat);
        Eigen::MatrixXd grad_P_mat = backprop_filter(sens_P_mat, rho_bar, beta);

        // 约束 2: 悬垂特征 Q
        Eigen::MatrixXd sens_Q_mat = Eigen::MatrixXd::Zero(rows, cols);
        double Q_val = compute_hanging_constraint(sens_Q_mat);
        Eigen::MatrixXd grad_Q_mat = backprop_filter(sens_Q_mat, rho_bar, beta);

        // 3. 体积 Vol
        Eigen::MatrixXd sens_V_mat = Eigen::MatrixXd::Zero(rows, cols);
        double V_val = compute_volume_constraint(sens_V_mat);
        Eigen::MatrixXd grad_V_mat = backprop_filter(sens_V_mat, apply_linear_filter(rho, cfg.r_min), beta);

        // 填入约束向量
        constraints(0) = P_val;
        constraints(1) = Q_val;
        constraints(2) = V_val;

        // 填入约束梯度矩阵 (Row 0 is P, Row 1 is Q)
        // 注意：这里需要确保 flatten 的顺序和 x 的顺序完全一致
        dfdx.row(0) = flatten(grad_P_mat);
        dfdx.row(1) = flatten(grad_Q_mat);
        dfdx.row(2) = flatten(grad_V_mat);

        mma.Update(x.data(), df0dx.data(), constraints.data(), dfdx.data(), xmin.data(), xmax.data());
        imposeSymmetry(x);

        // 注意: jdumas/mma 的 Update 函数对于 dgdx 通常期望一个转置的、或者特定顺序的指针
        // 标准 MMA 实现中 dgdx 是 m x n 的矩阵。
        // 为了安全起见，我们将 Eigen Matrix 转为 std::vector 或直接用 data()，
        // 但 Eigen 默认是 Column-Major，而 dgdx 通常期望 Constraint-Major (Row-Major 展平)。
        // 所以我们需要手动转置一下或者在 flatten 时注意。
        // 最稳妥的方式：构建一个展平的 vector，顺序为 [dP/dx1, dP/dx2... , dQ/dx1, dQ/dx2...]
        //std::vector<double> dgdx_vec(m * n);
        //Eigen::VectorXd flat_P = flatten(grad_P_mat);
        //Eigen::VectorXd flat_Q = flatten(grad_Q_mat);

        //// 拼接梯度
        //std::copy(flat_P.data(), flat_P.data() + n, dgdx_vec.begin());
        //std::copy(flat_Q.data(), flat_Q.data() + n, dgdx_vec.begin() + n);

        //// --- D. 调用 MMA 更新 ---
        //// Update(x, df0dx, f, dfdx, xmin, xmax)
        //// 注意：x 在这里既是输入也是输出
        //mma.Update(x.data(), df0dx.data(), constraints.data(), dgdx_vec.data(), xmin.data(), xmax.data());

        // --- E. 收敛检查与输出 ---
        double change = (x - xold).cwiseAbs().maxCoeff();
        xold = x;

        if (iter % 10 == 0) {
            // 可选：中间保存图片
            std::cout << "Iter: " << iter
                << " | J: " << J_val
                << " | P: " << P_val
                << " | Q: " << Q_val
                << " | Vol: " << rho_tilde.mean()
                << " | Beta: " << beta
                << " | Change: " << change << std::endl;
            unflatten(x, rho); // 确保 rho 是最新的
            update_rho_tilde(); // 辅助函数，只为了更新 rho_tilde 用于绘图
            printTopology();
            //save_result("debug_" + std::to_string(iter) + ".png");
        }

        // --- F. Continuation Scheme (Beta 更新) ---
        // 当 MMA 找到当前 Beta 下的解（变化率很小，且约束满足）时，增大 Beta
        if (change < 0.02 && P_val < 0.05 && Q_val < 0.05 && iter > 20) {
            if (beta < cfg.beta_max) {
                beta = std::min(cfg.beta_max, beta * cfg.beta_increase);
                std::cout << ">>> Increasing Beta to " << beta << std::endl;
                // 增大 Beta 后，梯度会突变，通常不需要重置 MMA，但可能会震荡一下
            }
            else {
                std::cout << ">>> Converged!" << std::endl;
                break;
            }
        }
        //if (iter % 50 == 0) {
        //    if (beta < cfg.beta_max)
        //        beta = std::min(cfg.beta_max, beta * cfg.beta_increase);
        //    //std::cout << "Updating Beta to: " << beta << std::endl;
        //    //if (P_val > 0.001 || Q_val > 0.001) {
        //    //    current_penalty = std::min(cfg.max_penalty, current_penalty * cfg.penalty_increase);
        //    //}
        //    //printTopology();
        //    // 输出状态

        //}

    }

    // 保存最终结果
    unflatten(x, rho);
    update_rho_tilde();
    save_result("output_mma_final.png");
}

/*void TopologyOptimizer::solve() {
    beta = cfg.beta_start;
    double current_penalty = cfg.penalty_weight; // 动态罚权重
    printTopology();
    Eigen::MatrixXd rho_prev = rho;
    double change_tolerance = 0.0002; // 收敛容差 (0.5% 的变化)

    for (int iter = 0; iter < cfg.max_iter; ++iter) {
        // 1. 过滤与投影 (Filter Chain)
        Eigen::MatrixXd rho_bar = apply_linear_filter(rho, cfg.r_min);
        apply_projection(rho_bar, beta); // 更新 rho_tilde

        // 2. 计算目标函数灵敏度 (几何逼近误差) [cite: 556]
        // J = 0.5 * || rho_tilde - target ||^2
        // dJ/d_rho_tilde = rho_tilde - target
        Eigen::MatrixXd sens_J = 1.0 * (rho_tilde - target_shape);


        // 3. 计算悬垂角约束及灵敏度 (Eq. 34, 45, 52)
        Eigen::MatrixXd sens_P = Eigen::MatrixXd::Zero(rows, cols);
        double P_val = compute_overhang_constraint(sens_P);

        // 计算悬垂特征约束及灵敏度 (Eq. 35, 46, 54)
        Eigen::MatrixXd sens_Q = Eigen::MatrixXd::Zero(rows, cols);
        double Q_val = compute_hanging_constraint(sens_Q);

        // 4. 自适应更新罚权重
    // 如果约束严重违背 (比如 > 0.01)，则增大惩罚力度
        if (iter % 10 == 0) {
            if (P_val > 0.001 || Q_val > 0.001) {
                current_penalty = std::min(cfg.max_penalty, current_penalty * cfg.penalty_increase);
            }
            double diff = (rho_tilde - target_shape).norm();
            std::cout << "Iter: " << iter << " | Diff: " << diff
                << " | P: " << P_val << " | Q: " << Q_val
                << " | Penalty: " << current_penalty  // 打印当前罚权重，观察是否在变大
                << " | Beta: " << beta << std::endl;
            printTopology();
        }

        // 5. 链式法则回传灵敏度 (Backprop through Filter)
        // dObj/drho = dObj/d_rho_tilde * d_rho_tilde/d_rho_bar * d_rho_bar/drho
        Eigen::MatrixXd total_grad = backprop_filter(sens_J, rho_bar, beta)
            + current_penalty * backprop_filter(sens_P, rho_bar, beta)
            + current_penalty * backprop_filter(sens_Q, rho_bar, beta);
        //Eigen::MatrixXd total_grad = backprop_filter(sens_J, rho_bar, beta)
        //    + 100 * backprop_filter(sens_P, rho_bar, beta)
        //    + 10 * backprop_filter(sens_Q, rho_bar, beta);

        double learning_rate = 0.05;
        double current_lr = learning_rate * std::exp(-0.002 * iter);

        // 设置一个下限，防止完全不动
        if (current_lr < 0.002) current_lr = 0.002;

        // 6. 更新设计变量 (简单的梯度下降)
        update_design_variables(total_grad, current_lr);
        imposeSymmetry();

        double max_change = (rho - rho_prev).cwiseAbs().maxCoeff();
        rho_prev = rho; // 更新旧值

        // 终止条件: 变化很小 且 约束基本满足
        // P_val 和 Q_val 是之前计算出的约束值
        if (max_change < change_tolerance && iter > 1000) {
            std::cout << ">>> CONVERGED at iter " << iter << "! Max change < " << change_tolerance << std::endl;
            break;
        }

        // 7. 更新 Beta (Continuation scheme)
        if (iter % 50 == 0) {
            if (beta < cfg.beta_max)
                beta = std::min(cfg.beta_max, beta * cfg.beta_increase);
            //std::cout << "Updating Beta to: " << beta << std::endl;
            //if (P_val > 0.001 || Q_val > 0.001) {
            //    current_penalty = std::min(cfg.max_penalty, current_penalty * cfg.penalty_increase);
            //}
            //printTopology();
            // 输出状态

        }



    }

    save_result("output_optimized.png");
}*/

void TopologyOptimizer::update_rho_tilde() {
    Eigen::MatrixXd rho_bar = apply_linear_filter(rho, cfg.r_min);
    apply_projection(rho_bar, beta, cfg.eta);
}

void TopologyOptimizer::printTopology() {
    std::cout << "\n--- Current Topology (Beta=" << beta << ") ---\n";
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double val = rho_tilde(i, j);
            //double val = rho_tilde(i, j);
            if (val > 0.75) std::cout << "#";      // 实体
            else if (val > 0.5) std::cout << "x"; // 灰色过渡带
            else if (val > 0.25) std::cout << "+"; // 灰色过渡带
            else std::cout << ".";                 // 空洞
        }
        std::cout << "\n";
    }
    std::cout << "-----------------------------------\n";
}


    // ==========================================
    // 3. 过滤器实现
    // ==========================================

Eigen::MatrixXd TopologyOptimizer::apply_linear_filter(const Eigen::MatrixXd& x, double r) {
    Eigen::MatrixXd x_new = Eigen::MatrixXd::Zero(rows, cols);
    Eigen::MatrixXd weight_sum = Eigen::MatrixXd::Zero(rows, cols);
    int r_int = std::ceil(r);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (int dy = -r_int; dy <= r_int; ++dy) {
                for (int dx = -r_int; dx <= r_int; ++dx) {
                    int ni = i + dy;
                    int nj = j + dx;
                    if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
                        double dist = std::sqrt(dx * dx + dy * dy);
                        if (dist <= r) {
                            double w = r - dist;
                            x_new(i, j) += w * x(ni, nj);
                            weight_sum(i, j) += w;
                        }
                    }
                }
            }
            x_new(i, j) /= (weight_sum(i, j) + 1e-10);
        }
    }
    return x_new;
}
void TopologyOptimizer::imposeSymmetry(Eigen::VectorXd & rho_vec) {
    // 假设设计域关于垂直中心线对称 (左右对称)
    //for (int i = 0; i < rows; ++i) {
    //    for (int j = 0; j < cols / 2; ++j) { // 只遍历左半边

    //        // 获取左侧和对称的右侧的值
    //        // Eigen 使用 (row, col) 访问
    //        double val_left = rho(i, j);
    //        double val_right = rho(i, cols - 1 - j);

    //        // 取平均值
    //        double avg_val = (val_left + val_right) / 2.0;

    //        // 强行赋值给两边
    //        rho(i, j) = avg_val;
    //        rho(i, cols - 1 - j) = avg_val;
    //    }
    //}
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols / 2; ++j) { // 只遍历左半边

            // 1. 手动计算一维索引 (假设之前是按行优先 flatten 的)
            // 左侧点 (i, j)
            int idx_left = i * cols + j;

            // 对称的右侧点 (i, cols - 1 - j)
            int idx_right = i * cols + (cols - 1 - j);

            // 2. 获取值 (使用 () 或 [] 访问一维向量)
            double val_left = rho_vec(idx_left);
            double val_right = rho_vec(idx_right);

            // 3. 取平均值
            double avg_val = (val_left + val_right) / 2.0;

            // 4. 强行赋值给两边
            rho_vec(idx_left) = avg_val;
            rho_vec(idx_right) = avg_val;
        }
    }
}

void TopologyOptimizer::apply_projection(const Eigen::MatrixXd& x_bar, double beta, double pro_thres) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            rho_tilde(i, j) = proj_heaviside(x_bar(i, j), beta, pro_thres);
        }
    }
}

Eigen::MatrixXd TopologyOptimizer::backprop_filter(const Eigen::MatrixXd& sens_phys, const Eigen::MatrixXd& rho_bar, double beta) {
    // 第一步：通过 Heaviside 投影的反向传播
    Eigen::MatrixXd sens_bar = Eigen::MatrixXd::Zero(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            sens_bar(i, j) = sens_phys(i, j) * d_proj_heaviside(rho_bar(i, j), beta, cfg.eta);
        }
    }
    // 第二步：通过线性滤波的反向传播 (线性滤波是对称算子，转置等于自身)
    return apply_linear_filter(sens_bar, cfg.r_min);
}

// ==========================================
// 4. 悬垂角约束逻辑 (The "Four Elements Pattern")
// ==========================================

// 获取"四元素模式"的全局索引 [cite: 281]
// 模式说明：以 (i,j) 为基准，Left-downward 模式涉及:
// 1(Self), 2(Left), 3(Left-Down), 4(Down)
// 注意：这里的方向依赖于你的网格定义。假设 Row 增加是向下。
// 左边是 Col-1， 下边是 Row+1。
void TopologyOptimizer::get_local_indices(int r, int c, std::vector<std::pair<int, int>>& idx_l, std::vector<std::pair<int, int>>& idx_r) {
    idx_l.clear(); idx_r.clear();

    // Left-Downward Pattern (Eq. 21 context)
    // 1: (0,0), 2: (-1,0) x轴向左, 3: (-1,-1), 4: (0,-1) y轴向下
    // 对应图像坐标：
    // 1: (r, c)
    // 2: (r, c-1)
    // 3: (r+1, c-1)
    // 4: (r+1, c)
    idx_l = { {r, c}, {r, c - 1}, {r + 1, c - 1}, {r + 1, c} };

    // Right-Downward Pattern
    // 1: (r, c)
    // 4: (r+1, c)
    // 5: (r+1, c+1)
    // 6: (r, c+1)
    // 对应论文图 7a 右侧
    idx_r = { {r, c}, {r + 1, c}, {r + 1, c + 1}, {r, c + 1} };
}

double TopologyOptimizer::compute_overhang_constraint(Eigen::MatrixXd& sensitivity) {
    double total_P = 0.0;
    double total_vol = rows * cols; // 简化体积计算
    double cos_crit = std::cos(cfg.alpha_crit_deg * M_PI / 180.0);
    double eps1 = 0.05; // [cite: 642]

    // 遍历每个元素计算约束 P (Eq. 45)
    for (int i = 0; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {

            std::vector<std::pair<int, int>> idx_l, idx_r;
            get_local_indices(i, j, idx_l, idx_r);

            // 提取局部密度
            Eigen::Vector4d rho_l, rho_r;
            for (int k = 0; k < 4; ++k) rho_l(k) = rho_tilde(idx_l[k].first, idx_l[k].second);
            for (int k = 0; k < 4; ++k) rho_r(k) = rho_tilde(idx_r[k].first, idx_r[k].second);

            // 计算梯度 [cite: 263, 286]
            // B_l 的前两行对应 dx, dy (梯度)
            Eigen::Vector2d grad_l = (B_l * rho_l).head<2>();
            Eigen::Vector2d grad_r = (B_r * rho_r).head<2>();

            double norm_l = grad_l.norm() + 1e-8; // 避免除零
            double norm_r = grad_r.norm() + 1e-8;

            // 计算 t 值 (Eq. 34)
            // t = grad^T * n - cos(alpha) * ||grad||
            double t_l = grad_l.dot(cfg.build_dir) - cos_crit * norm_l;
            double t_r = grad_r.dot(cfg.build_dir) - cos_crit * norm_r;

            // Sigmoid 激活 [cite: 504]
            double h_l = sigmoid(t_l - eps1);
            double h_r = sigmoid(t_r - eps1);
            double lambda = h_l * h_r; // 只有左右都违反时才惩罚

            // 累加约束值 P_i = rho * lambda * v
            double rho_i = rho_tilde(i, j);
            total_P += rho_i * lambda;

            // --- 计算灵敏度 (Eq. 52) ---
            // dP/drho_k = (Term1 + Term2)
            // 我们在这里使用"累加法"：计算当前元素 i 对其所有邻居 k 的梯度贡献

            // 1. 对自身的直接导数: d(rho_i * lambda)/drho_i = lambda + rho_i * dlambda/drho_i
            // 但 lambda 还依赖于 neighbors。

            double dlambda_dtl = d_sigmoid(t_l - eps1) * h_r;
            double dlambda_dtr = h_l * d_sigmoid(t_r - eps1);

            // 计算 dt/dgrad
            // dt/dgrad = n - cos * (grad / norm)
            Eigen::Vector2d dt_dgrad_l = cfg.build_dir - cos_crit * (grad_l / norm_l);
            Eigen::Vector2d dt_dgrad_r = cfg.build_dir - cos_crit * (grad_r / norm_r);

            // 贡献回传给 Left pattern 的 4 个邻居
            for (int k = 0; k < 4; ++k) {
                int ni = idx_l[k].first;
                int nj = idx_l[k].second;

                // dgrad_l / drho_k 是 B_l 的第 k 列 (前两行)
                Eigen::Vector2d dgrad_drho = B_l.col(k).head<2>();
                double dt_drho = dt_dgrad_l.dot(dgrad_drho);

                // 链式法则贡献
                // P_i = rho_i * lambda
                // dP_i / drho_k (via lambda) = rho_i * (dlambda/dtl * dtl/drho_k)
                sensitivity(ni, nj) += rho_i * dlambda_dtl * dt_drho / total_vol;
            }

            // 贡献回传给 Right pattern 的 4 个邻居
            for (int k = 0; k < 4; ++k) {
                int ni = idx_r[k].first;
                int nj = idx_r[k].second;
                Eigen::Vector2d dgrad_drho = B_r.col(k).head<2>();
                double dt_drho = dt_dgrad_r.dot(dgrad_drho);
                sensitivity(ni, nj) += rho_i * dlambda_dtr * dt_drho / total_vol;
            }

            // 别忘了 P 对 rho_i 自身的直接项 (Product rule)
            sensitivity(i, j) += lambda / total_vol;
        }
    }
    return total_P / total_vol;
}

// ==========================================
// 5. 悬垂特征约束逻辑 (Hanging Feature)
// ==========================================

double TopologyOptimizer::compute_hanging_constraint(Eigen::MatrixXd& sensitivity) {
    double total_Q = 0.0;
    double total_vol = rows * cols;
    int m = cfg.hang_m;
    double eps2 = 0.02; // [cite: 642]

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // 水平邻域平均密度 [cite: 421]
            // Left side check
            double sum_left = 0;
            int count_left = 0;
            // 注意：论文定义的模式是 2,4,6... 偶数位，或者简化为连续 m 个
            // 这里我们简化为连续 m 个邻居
            for (int k = 1; k <= m; ++k) {
                if (j - k >= 0) {
                    sum_left += rho_tilde(i, j - k);
                    count_left++;
                }
            }
            double avg_left = (count_left > 0) ? sum_left / m : 0; // 论文公式分母固定为 m

            // Right side check
            double sum_right = 0;
            int count_right = 0;
            for (int k = 1; k <= m; ++k) {
                if (j + k < cols) {
                    sum_right += rho_tilde(i, j + k);
                    count_right++;
                }
            }
            double avg_right = (count_right > 0) ? sum_right / m : 0;

            // Tau calculation [cite: 35, 36]
            double tau_l = rho_tilde(i, j) - avg_left - eps2;
            double tau_r = rho_tilde(i, j) - avg_right - eps2;

            // Sigmoid aggregation [cite: 505]
            // 只有当左右两侧都空（tau > 0）时，gamma 才为 1，即孤立点
            double h_l = sigmoid(tau_l);
            double h_r = sigmoid(tau_r);
            double gamma = h_l * h_r;

            double rho_i = rho_tilde(i, j);
            total_Q += rho_i * gamma;

            // --- 灵敏度计算 (Eq. 54) ---
            // Q_i = rho_i * gamma
            // 1. 直接项：dQi/drho_i = gamma + rho_i * (dgamma/drho_i)
            // dgamma/drho_i = dh_l * (1) * h_r + h_l * dh_r * (1)  (因为 dtau/drho_i = 1)

            double dgamma_dtaul = d_sigmoid(tau_l) * h_r;
            double dgamma_dtaur = h_l * d_sigmoid(tau_r);

            sensitivity(i, j) += (gamma + rho_i * (dgamma_dtaul + dgamma_dtaur)) / total_vol;

            // 2. 邻居项：dQi/drho_neighbor
            // Left neighbors (k from 1 to m)
            // dtau_l / drho_neigh = -1/m
            for (int k = 1; k <= m; ++k) {
                if (j - k >= 0) {
                    // 影响左侧邻居
                    sensitivity(i, j - k) += rho_i * dgamma_dtaul * (-1.0 / m) / total_vol;
                }
            }
            // Right neighbors
            for (int k = 1; k <= m; ++k) {
                if (j + k < cols) {
                    sensitivity(i, j + k) += rho_i * dgamma_dtaur * (-1.0 / m) / total_vol;
                }
            }
        }
    }
    return total_Q / total_vol;
}

double TopologyOptimizer::compute_volume_constraint(Eigen::MatrixXd& sensitivity) {
    double current_vol = rho_tilde.mean(); // mean = sum / N
    double max_vol = cfg.vol_frac; // 允许最大 50% 体积

    // g = V - V_max
    double val = current_vol - max_vol;

    // d(Vol)/d(rho_tilde) = 1/N
    sensitivity = Eigen::MatrixXd::Constant(rows, cols, 1.0 / (rows * cols));

    return val;
}
// ==========================================
// 6. 更新与IO
// ==========================================

//void TopologyOptimizer::update_design_variables(const Eigen::MatrixXd& grad, double learning_rate) {
//    // 简单的投影梯度下降
//    double damping = 0.5;
//
//    Eigen::MatrixXd rho_pro = rho - learning_rate * grad;
//    //rho = rho - learning_rate * grad;
//        
//    // 投影到 [0, 1]
//    for (int i = 0; i < rows; ++i) {
//        for (int j = 0; j < cols; ++j) {
//            // 简单的 Move limit
//            //double old_val = rho_tilde(i, j); // 近似
//            // Clamp
//            if (rho_pro(i, j) > 1.0) rho_pro(i, j) = 1.0;
//            if (rho_pro(i, j) < 0.0) rho_pro(i, j) = 0.0;
//        }
//    }
//    rho = damping * rho + (1.0 - damping) * rho_pro;
//
//}

void TopologyOptimizer::update_design_variables(const Eigen::MatrixXd& grad, double learning_rate) {
    double momentum = 0.9;      // 动量系数 (0.5 ~ 0.9)
    double move_limit = 0.05;

    // 速度更新公式: v = m * v - lr * grad
    velocity = momentum * velocity - learning_rate * grad;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double change = velocity(i, j);

            // 依然保留 Move Limit 作为最后一道防线
            if (change > move_limit) change = move_limit;
            if (change < -move_limit) change = -move_limit;

            double new_val = rho(i, j) + change;

            // 投影 [0, 1]
            if (new_val > 1.0) new_val = 1.0;
            if (new_val < 0.0) new_val = 0.0;

            // 反向修正速度（如果撞墙了，就把速度归零，防止堆积）
            if (new_val == 1.0 || new_val == 0.0) {
                velocity(i, j) = 0.0;
            }

            rho(i, j) = new_val;
        }
    }
}

void TopologyOptimizer::save_result(const std::string& filename) {
    cv::Mat out_img(rows, cols, CV_8UC1);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            out_img.at<uint8_t>(i, j) = static_cast<uint8_t>(rho_tilde(i, j) * 255.0);
        }
    }
    cv::imwrite(filename, out_img);
}

    