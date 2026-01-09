#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>

//#include "MMASolver.h"
// ==========================================
// 1. 配置与辅助结构
// ==========================================

#define M_PI 3.1415926

struct Config {
    int max_iter = 10000;
    double vol_frac = 0.5;      // 目标不需要体积约束，但保留参数习惯
    double r_min = 2.0;         // 滤波半径
    double alpha_crit_deg = 45.0; // 临界悬垂角
    int hang_m = 3;             // 悬垂特征水平检测长度
    double beta_start = 1.0;    // Heaviside投影参数
    double beta_max = 24.0;
    double move_limit = 0.1;    // 每次迭代最大变化量

    double penalty_weight = 10.0;    // 初始罚权重，从10改为100
    double penalty_increase = 1.1;    // 每次迭代如果没有满足约束，权重增加10%
    double max_penalty = 600.0;     // 权重上限，防止数值爆炸
    double eta = 0.3;
	double beta_increase = 1.2;
    // 构建方向：假设图像坐标系 (row, col)，Y向下，X向右。
    // 如果构建方向是"向上打印"，则法向量指向 (0, -1)
    Eigen::Vector2d build_dir = Eigen::Vector2d(0.0, -1.0);
};


// Sigmoid / Heaviside 函数及其导数
inline double sigmoid(double x, double k = 50.0) {
    return 1.0 / (1.0 + std::exp(-k * x));
}

inline double d_sigmoid(double x, double k = 50.0) {
    double s = sigmoid(x, k);
    return k * s * (1.0 - s);
}

// 投影函数 (Density Filter Heaviside)
inline double proj_heaviside(double rho, double beta, double eta = 0.5) {
    double num = std::tanh(beta * eta) + std::tanh(beta * (rho - eta));
    double den = std::tanh(beta * eta) + std::tanh(beta * (1.0 - eta));
    return num / den;
}

inline double d_proj_heaviside(double rho, double beta, double eta = 0.5) {
    double den = std::tanh(beta * eta) + std::tanh(beta * (1.0 - eta));
    double sec = 1.0 / std::cosh(beta * (rho - eta));
    return (beta * sec * sec) / den;
}


void initTShape();


// ==========================================
// 2. 拓扑优化核心类
// ==========================================


class TopologyOptimizer {
public:
    int rows, cols;
    double beta;
    Eigen::MatrixXd rho;          // 设计变量 (0-1)
    Eigen::MatrixXd rho_tilde;    // 物理密度 (滤波+投影后)
    Eigen::MatrixXd target_shape; // 初始目标形状 (用于几何逼近目标函数)
    Config cfg;


    // 固定的局部拟合矩阵 B_l (Ref: Eq. 25)
    Eigen::Matrix<double, 3, 4> B_l;
    Eigen::Matrix<double, 3, 4> B_r;

    TopologyOptimizer(const std::string& image_path, Config config);

    void solve();

    void printTopology();

    // ==========================================
    // 3. 过滤器实现
    // ==========================================

    Eigen::MatrixXd apply_linear_filter(const Eigen::MatrixXd& x, double r);

    void imposeSymmetry();

    void apply_projection(const Eigen::MatrixXd& x_bar, double beta);

    Eigen::MatrixXd backprop_filter(const Eigen::MatrixXd& sens_phys, const Eigen::MatrixXd& rho_bar, double beta);

    // ==========================================
    // 4. 悬垂角约束逻辑 (The "Four Elements Pattern")
    // ==========================================

    // 获取"四元素模式"的全局索引 [cite: 281]
    // 模式说明：以 (i,j) 为基准，Left-downward 模式涉及:
    // 1(Self), 2(Left), 3(Left-Down), 4(Down)
    // 注意：这里的方向依赖于你的网格定义。假设 Row 增加是向下。
    // 左边是 Col-1， 下边是 Row+1。
    void get_local_indices(int r, int c, std::vector<std::pair<int, int>>& idx_l, std::vector<std::pair<int, int>>& idx_r);

    double compute_overhang_constraint(Eigen::MatrixXd& sensitivity);

    // ==========================================
    // 5. 悬垂特征约束逻辑 (Hanging Feature)
    // ==========================================

    double compute_hanging_constraint(Eigen::MatrixXd& sensitivity);

    // ==========================================
    // 6. 更新与IO
    // ==========================================

    void update_design_variables(const Eigen::MatrixXd& grad, double learning_rate);

    void save_result(const std::string& filename);

};

