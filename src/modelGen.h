#pragma once
#include "Tool.h"
#include <queue>
#include <random>
#include <chrono>  // 添加时间测量

// 定义高斯核的数据结构
class GaussianKernel {

public:
    GaussianKernel() : center(Eigen::Vector3d::Zero()), sigma(0.1), amplitude(1.0) , center_value(0.0), on_surface(false){}
    /*GaussianKernel() {};*/
    GaussianKernel(Eigen::Vector3d cente_, double sigma_, double amplitude_, double center_value_ = 0.0);
	
    double gaussian_fun(const Eigen::Vector3d& p);
    bool is_on_surface() const;

public:
    Eigen::Vector3d center; // 核的中心位置
    double sigma;         // 核的大小/影响力范围 (高斯函数的标准差)
    double amplitude;
    double center_value;
    bool on_surface;

};

struct Edge {
    int from;
    int to;
    double length;
    double weight;
};

struct CompareEdge {
    bool operator()(const Edge& e1, const Edge& e2) const {
        return e1.weight > e2.weight;
    }
};

using Graph = std::vector<std::vector<Edge>>;

class ModelGenerator {
public:
    ModelGenerator() {};
    ModelGenerator(std::string input_file, int pores = PoresNum);

	void generateGaussianSDF();

    void sample_interior_points(std::vector<Eigen::Vector3d>& pore_centers, std::vector<double>& pore_sdfs, 
        std::vector<int>& inside_indices, int pores, std::mt19937& gen);

    void generate_gaussians(std::vector<Eigen::Vector3d> pore_centers, std::vector<double> pore_sdfs, std::mt19937& gen);

    double combinedSDF(Eigen::Vector3d& p, std::vector<GaussianKernel> G_kernels, double C);

    void show_model();

    std::vector<Edge>  pores_connection_mst(const std::vector<GaussianKernel>& gau);
    std::vector<int> all_leafs_mst(std::vector<Edge>& mst_tree);

    double generate_tube(const Eigen::Vector3d& p, const GaussianKernel& k1, const GaussianKernel& k2, double iso_level_C, double mid_radius_factor);    
    double generate_tube2( Eigen::Vector3d& p,  GaussianKernel& k1,  GaussianKernel& k2, double iso_level_C, double mid_radius_factor = 0.5);

    double calculate_edge_weight(GaussianKernel k1, GaussianKernel k2);

    std::vector<int> find_path_in_tree(int start_node_id, int end_node_id, std::vector<Edge> graph, int num_nodes);
	double length_graph_path(int p1, int p2);
    double length_path(int p1, int p2);
    int find_edge_by_nodes(int from_node, int to_node, const std::vector<Edge> edge_list);

    std::pair<double, double> calculate_each_path(const std::vector<int>& path);
    double calculate_score(std::vector<std::vector<int>>  Paths);

private:
    double m_currentPorosity = 0; 
    // 缓存当前可视化模型对应的SDF与网格点，以便独立后处理
	int pore_num = PoresNum;			   // 空洞数量
    int resolution = Resolution;               // 网格分辨率
    double isolevel = Isolevel;
    double gauss_iso = Gauss_level;

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
    double safe_distance = 0; //dart throwing 
    double amplitude_min = Amplitude_min;
    double amplitude_max = Amplitude_max;
    double sigma_min = Sigma_min;
    double sigma_max = Sigma_max;
    std::vector<GaussianKernel> Kernels;
    std::vector<int> surface_kernels;

    std::vector<Edge> Tube_edges;


	double smooth_t = SmoothT;         //平滑参数，值越大，平滑效果越小，趋近于普通并集
    int m_cachedRes = 0;                   // 分辨率
    double cut_face = 0.25;             //切割底座
    bool m_sdfValid = false;               // 缓存有效标志
    Eigen::VectorXd m_originalCachedSDF;   // 初始（未后处理）SDF备份
    bool m_hasPostProcessed = false;       // 是否已经做过后处理
};