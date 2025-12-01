#pragma once
#include "Tool.h"
#include <queue>
#include <random>
#include <chrono>  // 添加时间测量

#include <filesystem>
#include <igl/read_triangle_mesh.h>

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

    // 重载 > 运算符以实现最小堆
    bool operator>(const Edge& other) const {
		if (Direct_dis) return length > other.length;
        else return weight > other.weight;
    }
};



struct NodeDist {
    int id;
    double dist;
    // 重载 > 运算符以实现最小堆
    bool operator>(const NodeDist& other) const {
        return dist > other.dist;
    }
};

using AdjacencyList = std::vector<std::vector<int>>;

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

    std::vector<Edge>  pores_connection_mst(const std::vector<GaussianKernel>& gau, int Dmax = 7);
    std::vector<std::vector<int>> construct_adj_list(std::vector<Edge> edges_list, int kernel_num);
    std::vector<double> construct_dist_map(int p_index, AdjacencyList adj);
    std::vector<std::vector<int>> get_unused_edge_adj(AdjacencyList Adj_list, double dis_thres);

    std::vector<int> find_path_in_tree(int p1, int p2, int num_nodes, AdjacencyList adj);
    double length_graph_path(int p1, int p2, AdjacencyList adj);
    double length_path(int p1, int p2);
    int find_edge_by_nodes(int from_node, int to_node, const std::vector<Edge> edge_list);
    bool find_edge_in_path(Edge cand_edge, vector<int> path);
    std::vector<int>  find_specified_path(int p_index, int s1, int s2, AdjacencyList adj, bool show_debug = false); //经过点p_index的，两端点为s1s2的路径

    std::vector<int> all_leafs_mst(std::vector<Edge>& mst_tree);

    double generate_tube(const Eigen::Vector3d& p, const GaussianKernel& k1, const GaussianKernel& k2, double iso_level_C, double mid_radius_factor);    
    double generate_tube2( Eigen::Vector3d& p,  GaussianKernel& k1,  GaussianKernel& k2, double iso_level_C, double mid_radius_factor = 0.5);

    double calculate_edge_weight(GaussianKernel k1, GaussianKernel k2);


	double calculate_path_translucency(std::vector<int>& path, bool show_debug = false);  //path里存放的是kernel的索引
    double cal_kernel_translucency(int p_index, int& max_s1, int& max_s2, std::vector<int>& max_path, AdjacencyList adj, bool debug=false);
    double cal_total_translucency(std::vector<GaussianKernel> gau, std::vector<int> surface_ks, AdjacencyList adj);

    vector<int> check_inner_leafs(vector<int> leafs_index);
	
    int find_nearest_grid(Eigen::Vector3d point);
	double line_cross_surface(Eigen::Vector3d p1, Eigen::Vector3d p2, double thres, int sam_num);


    double calculate_score(std::vector<std::vector<int>>  Paths);

    //---------------optimize------------------
    vector<int> cal_edge_usage(std::vector<std::vector<int>> Paths);
    pair<double, double> add_edges(Edge cand_edge, AdjacencyList adj, std::vector<int>& max_path1, std::vector<int>& max_path2);
    bool replace_edges(int p_index, int replace_e, std::vector<Edge>& Tube_edges, AdjacencyList& adj, AdjacencyList& unused_adj);
    void optimize_mst(int opt_times_once, int edge_max, bool debug = false);

	//------------generate tubes----------------
    int generate_mst_tubes(int grid_num, int res, double iso, double gaus_iso, double smooth_t);

private:
    double m_currentPorosity = 0; 
    // 缓存当前可视化模型对应的SDF与网格点，以便独立后处理
	int pore_num = PoresNum;			   // 空洞数量
    int resolution = Resolution;               // 网格分辨率

    Eigen::MatrixXd V_ini; //初始网格顶点
    Eigen::MatrixXi F_ini; // 初始网格面片
    Eigen::VectorXd SDF_ini;           // 原始/已处理的SDF网格值
    Eigen::MatrixXd GV;            // 网格点坐标

    Eigen::MatrixXd V_out; //输出网格顶点
    Eigen::MatrixXi F_out; // 输出网格面片
    Eigen::VectorXd SDF_out;           // 输出的SDF网格值
    //Eigen::MatrixXd GV_out;            // 网格点坐标


	//高斯核参数范围
    double safe_distance = 0; //dart throwing 
    double amplitude_min = Amplitude_min;
    double amplitude_max = Amplitude_max;
    double sigma_min = Sigma_min;
    double sigma_max = Sigma_max;
    std::vector<GaussianKernel> Kernels;
    std::vector<int> surface_kernels;
    std::vector<std::vector<int>> Paths;
    vector<pair<int, int>> max_paths_kernel;  //每个kernel通透性最大的路径两端
    vector<double> kernel_translucency;

    std::vector<Edge> Tube_edges;
    std::vector<std::vector<int>> Adj_list;
    std::vector<std::vector<int>> Unused_adj_list;


    double finalPorosity = 0;
	double smooth_t = SmoothT;         //平滑参数，值越大，平滑效果越小，趋近于普通并集
    int m_cachedRes = 0;                   // 分辨率
    double cut_face = 0.25;             //切割底座
    bool m_sdfValid = false;               // 缓存有效标志
    Eigen::VectorXd m_originalCachedSDF;   // 初始（未后处理）SDF备份
    bool m_hasPostProcessed = false;       // 是否已经做过后处理
};