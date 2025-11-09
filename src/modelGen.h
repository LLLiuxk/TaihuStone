#pragma once
#include "Tool.h"
#include <random>
#include <chrono>  // 添加时间测量

// 定义高斯核的数据结构
class GaussianKernel {

public:
    GaussianKernel() {};
    GaussianKernel(Eigen::Vector3d cente_, double sigma_, double amplitude_);
	//GaussianKernel() : center(Eigen::Vector3d::Zero()), sigma(0.1), amplitude(1.0) {}
    double gaussian_fun(const Eigen::Vector3d& p);

public:
    Eigen::Vector3d center; // 核的中心位置
    double sigma;         // 核的大小/影响力范围 (高斯函数的标准差)
    double amplitude;

};



class ModelGenerator {
public:
    ModelGenerator() {};
    ModelGenerator(std::string input_file, int pores = PoresNum);

	void generateGaussianSDF();
    double combinedSDF(Eigen::Vector3d& p);

    void show_model();

private:
    double m_currentPorosity = 0; 
    // 缓存当前可视化模型对应的SDF与网格点，以便独立后处理
	int pore_num = PoresNum;			   // 空洞数量
    int resolution = Resolution;               // 网格分辨率
    double isolevel = Isolevel;
    double gauss_combined = Gauss_level;

    Eigen::MatrixXd V_ini; //初始网格顶点
    Eigen::MatrixXi F_ini; // 初始网格面片
    Eigen::VectorXd SDF_ini;           // 原始/已处理的SDF网格值
    Eigen::MatrixXd GV;            // 网格点坐标

    Eigen::MatrixXd V_out; //输出网格顶点
    Eigen::MatrixXi F_out; // 输出网格面片
    Eigen::VectorXd SDF_out;           // 输出的SDF网格值
    //Eigen::MatrixXd GV_out;            // 网格点坐标

    double finalPorosity = 0;

	//高斯核参数范围
    double safe_distance;
    double amplitude_min = Amplitude_min;
    double amplitude_max = Amplitude_max;
    double sigma_min = Sigma_min;
    double sigma_max = Sigma_max;
    std::vector<GaussianKernel> Kernels;

	double smooth_t = SmoothT;         //平滑参数，值越大，平滑效果越小，趋近于普通并集
    int m_cachedRes = 0;                   // 分辨率
    double cut_face = 0.25;             //切割底座
    bool m_sdfValid = false;               // 缓存有效标志
    Eigen::VectorXd m_originalCachedSDF;   // 初始（未后处理）SDF备份
    bool m_hasPostProcessed = false;       // 是否已经做过后处理
};