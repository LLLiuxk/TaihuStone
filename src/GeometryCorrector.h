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

extern double Match_weight;

extern double R_min;

extern double Vol_weight;

extern double Overhang_weight;

extern double Vol_frac;

extern double Max_beta;

struct Vec2 {
    double x, y;
};

class GeometryCorrector {
public:
    int nx, ny;                 // 网格分辨率 (宽, 高)
    double r_min;               // 滤波半径

    int nx_logical, ny_logical; // 逻辑分辨率 (设计意图的尺寸)
    int padding;

    std::vector<double> x;          // 原始设计变量 (Design Variables)
    std::vector<double> rho;       // 设计变量 (当前形状)
    std::vector<double> rho_init;  // 初始输入形状 (目标形状)
    std::vector<double> rho_tilde;  // 滤波后的密度 (Filtered Density)
    std::vector<double> sensitivity; // 目标函数的梯度 (dJ/drho)

    double beta; // 投影锐度参数 (逐渐增大以获得清晰边界)
    double eta;  // 投影阈值 (通常为 0.5)



public:
    GeometryCorrector(int logical_w, int logical_h, int pad, double filter_radius);

    // 初始化：创建一个简单的测试形状（例如一个悬空的圆，不可打印）
    void initializeTestShape();

    void initTShape();

    // 辅助函数：设置密度（带边界检查）
    void setDensity(std::vector<double>& field, int x, int y, double val);
    // 辅助函数：获取密度
    double getDensity(const std::vector<double>& field, int x, int y);

    // ---------------------------------------------------------
    // 核心逻辑 1: 计算目标函数 (几何偏差) 及其导数
    // J = Sum( (rho_i - rho_init_i)^2 )
    // dJ/drho = 2 * (rho_i - rho_init_i)
    // ---------------------------------------------------------
    double computeObjectiveAndGradient();

    // ---------------------------------------------------------
    // 核心逻辑 2: 密度滤波 (Density Filter)
    // 作用：平滑设计变量，引入长度尺度，避免棋盘格。
    // 公式：rho_tilde_i = Sum(w_ij * x_j) / Sum(w_ij)
    // ---------------------------------------------------------
    void applyFilter();
    // ---------------------------------------------------------
    // 核心逻辑 3: Heaviside 投影 (Projection)
    // 作用：将灰度推向 0 或 1。
    // 公式：使用 tanh 近似 Heaviside 函数
    // ---------------------------------------------------------
    void applyProjection();

    // ---------------------------------------------------------
    // 核心逻辑 4: 链式法则求导 (Chain Rule)
    // 目标：计算 dJ/dx (对设计变量的导数)，用于梯度下降。 
    // 路径：dJ/dx = (dJ/drho) * (drho/drho_tilde) * (drho_tilde/dx)
    // ---------------------------------------------------------
    void applyChainRule();
    // ---------------------------------------------------------
   // 核心逻辑 5: 悬垂约束 (Overhang Constraint)
   // 这是一个简化的几何约束：惩罚那些“上方有材料但下方无支撑”的单元。
   // ---------------------------------------------------------
    void applyOverhangConstraint(double penalty_weight);

    void resetBoundarySensitivities();

    void printTopology();

    // 获取设计变量引用，用于优化器更新
    std::vector<double>& getDesignVariables() { return x; }
    // 获取灵敏度引用
    const std::vector<double>& getSensitivities() { return sensitivity; }
    // 更新 Beta (用于延拓法)
    void increaseBeta() { if (beta < Max_beta) beta *= 2.0; }

    // ---------------------------------------------------------
    // 工具：导出为 PPM 图片 (无需库即可查看)
    // ---------------------------------------------------------
    void exportPPM(const std::string& filename, const std::vector<double>& field);

    void imposeSymmetry();
};

void gradientDescentUpdate(std::vector<double>& x,
    const std::vector<double>& dc,
    double learning_rate);

void optimalityCriteriaUpdate(std::vector<double>& x,
    const std::vector<double>& dc,
    double target_vol_frac);