#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "MMASolver.h"

struct Config3D {
    int max_iter = 500;        // 3D 计算量大，通常迭代次数少一些
    double vol_frac = 0.5;     // 体积约束
    double r_min = 2.0;        // 3D 滤波半径 (单位：体素)
    double alpha_crit_deg = 45.0;

    // 3D 悬挂检测通常看垂直投影或局部支撑
    // 这里 hang_m 代表水平邻域检测半径
    int hang_m = 2;

    double beta_start = 1.0;
    double beta_max = 64.0;    // 3D 建议不要设太高，梯度容易断
    double beta_increase = 1.5;

    double cost_remove = 1000.0;
    double cost_add = 0.001;

    double eta = 0.5;          // 投影阈值

    // 构建方向通常是 Z 轴 (0,0,1)
    Eigen::Vector3d build_dir = Eigen::Vector3d(0.0, 0.0, 1.0);
};

class TopologyOptimizer3D {
public:
    // 维度信息
    int nx, ny, nz; // x:宽, y:深, z:高
    int n_vars;     // nx * ny * nz

    Config3D cfg;
    double beta;

    // 核心数据场 (Flattened 1D Vectors)
    // 使用 Eigen::VectorXd 方便与 MMA 交互
    Eigen::VectorXd rho;          // 设计变量 x
    Eigen::VectorXd rho_tilde;    // 物理密度 (Projected)
    Eigen::VectorXd target_shape; // 目标形状 (Reference)

    // 辅助场
    Eigen::VectorXd sensitivity_map; // 用于可视化灵敏度

    // 构造函数：传入二进制文件路径
    TopologyOptimizer3D(const std::string& vol_path, Config3D config);

    // 主求解循环
    void solve();

    // ==========================================
    // 核心工具函数
    // ==========================================

    // 坐标索引映射: (x, y, z) -> idx
    inline int idx(int x, int y, int z) const {
        // 边界检查 (Debug模式下可用，Release建议去掉以提速)
        // if(x<0 || x>=nx || y<0 || y>=ny || z<0 || z>=nz) return -1;
        return z * (nx * ny) + y * nx + x;
    }

    // 逆映射 (仅用于调试)
    inline void idx2coord(int id, int& x, int& y, int& z) const {
        z = id / (nx * ny);
        int rem = id % (nx * ny);
        y = rem / nx;
        x = rem % nx;
    }

    // IO 相关
    void loadVolumeData(const std::string& path);
    void saveVTI(const std::string& filename); // 保存为 ParaView 格式

    // 物理场更新
    void update_physics_field();
    void imposeSymmetry(Eigen::VectorXd& x); // 3D 对称性

    // 滤波与投影
    // 3D 线性滤波 (显式卷积实现，虽然慢但简单)
    Eigen::VectorXd apply_linear_filter(const Eigen::VectorXd& x);

    // 滤波的反向传播
    Eigen::VectorXd backprop_filter(const Eigen::VectorXd& sens_phys, const Eigen::VectorXd& rho_bar);

    // 辅助数学函数
    double proj_heaviside(double val, double b, double e) {
        double num = std::tanh(b * e) + std::tanh(b * (val - e));
        double den = std::tanh(b * e) + std::tanh(b * (1.0 - e));
        return num / den;
    }

    double d_proj_heaviside(double val, double b, double e) {
        double den = std::tanh(b * e) + std::tanh(b * (1.0 - e));
        double t = std::tanh(b * (val - e));
        return b * (1.0 - t * t) / den;
    }

    // ==========================================
    // 待实现的约束函数 (将在下一部分详细展开)
    // ==========================================
    double compute_objective(Eigen::VectorXd& sens);
    double compute_volume_constraint(Eigen::VectorXd& sens);
    double compute_overhang_constraint(Eigen::VectorXd& sens); // 3D版
    double compute_hanging_constraint(Eigen::VectorXd& sens);  // 3D版

private:
    // 预计算的邻居列表，用于加速 3D 滤波
    // filter_neighbors[i] 存储了体素 i 周围所有在 r_min 内的邻居索引和权重
    struct Neighbor { int idx; double weight; };
    std::vector<std::vector<Neighbor>> filter_neighborhoods;
    std::vector<double> filter_weight_sums;

    void precompute_filter(); // 初始化时调用
};
