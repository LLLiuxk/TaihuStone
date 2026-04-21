#pragma once

#include "Tool.h"
#include <queue>
#include <random>
#include <chrono>  // 添加时间测量
#include <direct.h>
#include <filesystem>
#include <omp.h>
#include <igl/read_triangle_mesh.h>


// 定义高斯核的数据结构
//class GaussianKernel {
//
//public:
//    GaussianKernel() : center(Eigen::Vector3d::Zero()), sigma(0.1), amplitude(1.0) , center_value(0.0), on_surface(false){}
//    /*GaussianKernel() {};*/
//    GaussianKernel(Eigen::Vector3d cente_, double sigma_, double amplitude_, double center_value_ = 0.0);
//	
//    double gaussian_fun(const Eigen::Vector3d& p);
//    bool is_on_surface() const;
//
//public:
//    Eigen::Vector3d center; // 核的中心位置
//    double sigma;         // 核的大小/影响力范围 (高斯函数的标准差), sigma越大值越大，等值面的圆越大
//    double amplitude;
//    double center_value;
//    bool on_surface;
//
//};
struct ParaSensitivityStats {
    double sum_T_angle = 0.0;
    double sum_T_length = 0.0;
    double sum_T_location = 0.0;
    double sum_S_horiz = 0.0;
    double sum_weighted_vp = 0.0;
    int valid_path_num = 0;
};

struct KernelRandomSpec {
    double amp01 = 0.0;
    double sigma01 = 0.0;
    double orient_u01 = 0.0;
    double orient_v01 = 0.0;
    double axis01 = 0.0;
};

struct RenderParamSnapshot {
    std::string tag = "param";
    double amplitude_min = 1.0;
    double amplitude_max = 1.0;
    double sigma_min = 0.015;
    double sigma_max = 0.033;
    double w4sig_max = 3.0;
    double w4sig_min = 25.0;
    double axis_max_ratio = 1.7;
    double gauss_level = 0.5;
    double smooth_t = 15.0;
    double tube_radius_factor = 0.8;
};

class GaussianKernel {

public:

    GaussianKernel(
        const Eigen::Vector3d& center_,
        double sigma_,
        double sigma_parallel_,              // 沿主轴（建造方向）
        double sigma_perp_,                  // 垂直主轴
        const Eigen::Matrix3d& rotation_,    // 主轴方向
        double amplitude_,
        double center_value_ = 0.0
    );

    GaussianKernel(
        const Eigen::Vector3d& center_,
        double sigma_,              // 沿主轴（建造方向）
        double amplitude_,
        double center_value_ = 0.0
    );

    double gaussian_fun(const Eigen::Vector3d& p);
    bool is_on_surface();
    
    void rebuildInvSigma();

public:
    Eigen::Vector3d center;

    double sigma_parallel;   // σ∥
    double sigma_perp;       // σ⊥

    Eigen::Matrix3d R;       // 局部 → 世界 旋转矩阵
    Eigen::Matrix3d invSigma; // Σ^{-1}

    double sigma;   // ball
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

class PathQuery {
public:
    PathQuery(int num_nodes, const std::vector<std::vector<int>>& adj, int root);
    // 查询从 root 到任意目标 t 的路径
    std::vector<int> query_path(int t);

    // 一次 BFS 建立 parent 数组
    bool build_parent_tree();

public:
    int n;
    int root;
    bool isConnected;
    std::vector<std::vector<int>> adj_list;
    std::vector<int> parent;

  
};

using AdjacencyList = std::vector<std::vector<int>>;

class ModelGenerator {
public:
    ModelGenerator() {};
    ModelGenerator(std::string input_file, int pores = PoresNum);

    void model_porous_structure();
	void generateGaussianSDF();
    void supportFreeOpt();

    void sample_interior_points(std::vector<Eigen::Vector3d>& pore_centers, std::vector<double>& pore_sdfs, 
        std::vector<int>& inside_indices, int pores, std::mt19937& gen);
    void sample_interior_close(std::vector<Eigen::Vector3d>& pore_centers, std::vector<double>& pore_sdfs,
        std::vector<int>& inside_indices, int pores, std::mt19937& gen);
    void sample_regular(std::vector<Eigen::Vector3d>& pore_centers, std::vector<double>& pore_sdfs, std::vector<int>& inside_indices, int pores);

    void generate_gaussians(std::vector<Eigen::Vector3d> pore_centers, std::vector<double> pore_sdfs, std::mt19937& gen);
    //void generate_gaussians_iso(std::vector<Eigen::Vector3d> pore_centers, std::vector<double> pore_sdfs, std::mt19937& gen);

    double combinedSDF(Eigen::Vector3d& p, std::vector<GaussianKernel> G_kernels, double C);

    void show_model();

    std::vector<Edge>  pores_connection_mbdst(const std::vector<GaussianKernel>& gau, int Dmax = 7);
    std::vector<Edge>  pores_connection_mst(const std::vector<GaussianKernel>& gau);
    std::vector<std::vector<int>> construct_adj_list(std::vector<Edge> edges_list, int kernel_num);
    std::vector<std::vector<int>> construct_adj_list(std::vector<pair<int, int>> edges_list, int kernel_num);
    std::vector<std::vector<int>> get_unused_edge_adj(AdjacencyList Adj_list, double dis_thres);
    double cal_path_graph_length(std::vector<int> path_);

    std::vector<int> find_path_in_tree(int p1, int p2, int num_nodes, AdjacencyList adj);
    std::vector<int>  find_specified_path(int p_index, int s1, int s2, AdjacencyList adj); //经过点p_index的，两端点为s1s2的路径
    double length_graph_path(int p1, int p2, AdjacencyList adj);
    double length_path(int p1, int p2);
    int find_edge_by_nodes(int from_node, int to_node, const std::vector<Edge> edge_list);
    bool find_edge_in_path(Edge cand_edge, vector<int> path);

    std::vector<int> all_leafs_mst(std::vector<Edge>& mst_tree);

    double calculate_edge_weight(GaussianKernel k1, GaussianKernel k2);

	
	double calculate_path_translucency(std::vector<int>& path, bool show_debug = false);  //path里存放的是kernel的索引
    double calculate_path_translucency2(std::vector<int>& path, bool show_debug);  //old version
    double cal_kernel_translucency(int p_index, int& max_s1, int& max_s2, std::vector<int>& max_path, AdjacencyList adj, bool debug=false);
    double cal_total_translucency(std::vector<GaussianKernel> gau, AdjacencyList adj);
    ParaSensitivityStats cal_para_sensitivity(bool show_debug = false);


    vector<int> check_inner_leafs(vector<int> leafs_index);
	
    int find_nearest_grid(Eigen::Vector3d point);
	double line_cross_surface(Eigen::Vector3d p1, Eigen::Vector3d p2, double thres, int sam_num = 3);

    //---------------optimize------------------
    vector<int> cal_edge_usage(std::vector<std::vector<int>> Paths, bool show_debug = true);
    pair<double, double> add_edges(Edge cand_edge, AdjacencyList adj, std::vector<int>& max_path1, std::vector<int>& max_path2, bool debug = false);
    bool replace_edges(int p_index, int replace_e, std::vector<Edge>& Tube_edges, AdjacencyList& adj, AdjacencyList& unused_adj);
    void optimize_mst(int opt_times_once, int edge_max, vector<int>& rep_vec, bool debug = false);
	void optimize_mst2(int itea_max_times, int max_edge, bool iter_add = false, bool debug = false); //max_edge = 0代表最大边数递增，！=0代表固定最大边数

	//------------generate tubes----------------
    double generate_tube(const Eigen::Vector3d& p, const GaussianKernel& k1, const GaussianKernel& k2, double iso_level_C, double mid_radius_factor);
    double generate_tube2(Eigen::Vector3d& p, GaussianKernel& k1, GaussianKernel& k2, double iso_level_C, double mid_radius_factor = 0.5);
    double generate_tube3(Eigen::Vector3d& p, int k1_index, int k2_index, vector<Eigen::Matrix3d>& S_matrixs,
        Eigen::Matrix2d W_perp, double iso_level_C, double mid_radius_factor = 0.5);
    int generate_mbdst_tubes(std::vector<pair<int, int>> edge_con, int grid_num, int res, double iso, double gaus_iso, double smooth_t);
    void calculate_tube_matrixs(vector<Eigen::Matrix3d>& S_matrixs, std::vector<std::vector<Eigen::Matrix2d>>& W_perp_matrixs);
    void compare_msc(Eigen::VectorXd SDF_gaussian, int res, int grid_num, double smooth_t);

    void test_item();

	void optimize_model_py(std::string& filename, std::string& outfilename);

    int resolveOverlaps3D(
        std::vector<GaussianKernel>& kernels,
        int maxIters = 50,
        double isoValue = 0.5,
        double minSigmaPerp = 1e-4,
        double minSigmaParallel = 1e-4,
        double tol = 1e-5,                 // 投影分离容差（单位同坐标）
        double margin = 1e-4,               // 额外分离裕量（可设成 1e-4 等）
        double minScalePerUpdate = 0.90,   // 每次更新最多缩到多少（避免骤缩/抖动）
        bool verbose = false);

    void buildCandidateAxes(const GaussianKernel& a, const GaussianKernel& b, const Eigen::Vector3d& delta, std::vector<Eigen::Vector3d>& axesOut);

    double supportHalfLengthOnAxis(const GaussianKernel& k, const Eigen::Vector3d& n_unit, double chi);

    void removeUnusedNodes(std::vector<GaussianKernel>& k, std::vector<int>& surface_kernels, std::vector<std::pair<int, int>>& edges, std::vector<std::vector<int>>& adj);

    RenderParamSnapshot captureCurrentRenderParams(const std::string& tag = "param1") const;
    RenderParamSnapshot captureSecondRenderParams(const std::string& tag = "param2") const;
    void applyRenderParams(const RenderParamSnapshot& params);
    void build_kernels_from_specs(
        const std::vector<Eigen::Vector3d>& pore_centers,
        const std::vector<double>& pore_sdfs,
        const std::vector<KernelRandomSpec>& random_specs);
    void render_fixed_skeleton_variant(
        const std::vector<Eigen::Vector3d>& pore_centers,
        const std::vector<double>& pore_sdfs,
        const std::vector<KernelRandomSpec>& random_specs,
        const std::vector<std::pair<int, int>>& fixed_edges,
        const RenderParamSnapshot& params,
        const std::string& output_dir,
        bool rebuild_kernels = true);
    void render_fixed_skeleton_dual_params(
        const std::vector<Eigen::Vector3d>& pore_centers,
        const std::vector<double>& pore_sdfs,
        const std::vector<KernelRandomSpec>& random_specs,
        const std::vector<std::pair<int, int>>& fixed_edges);

private:

	int pore_num = PoresNum;			   // 空洞数量
    int resolution = Resolution;               // 网格分辨率

    Eigen::MatrixXd V_ini; //初始网格顶点
    Eigen::MatrixXi F_ini; // 初始网格面片
    Eigen::VectorXd SDF_ini;           // 原始/已处理的SDF网格值
    Eigen::MatrixXd GV;            // 网格点坐标
    Eigen::MatrixXd V_t; //输出高斯场对应网格顶点
    Eigen::MatrixXi F_t; // 输出高斯场对应网格面片

    Eigen::MatrixXd V_out; //输出网格顶点
    Eigen::MatrixXi F_out; // 输出网格面片
    Eigen::VectorXd SDF_out;           // 输出的SDF网格值
    Eigen::VectorXd SDF_gau;           // 只有gaussian的SDF网格值
    //Eigen::MatrixXd GV_out;            // 网格点坐标
    VoxelGrid Grids;

    Eigen::Vector3d bb_min; 
    Eigen::Vector3d bb_max;

	//高斯核参数范围
    double safe_distance = 0; //dart throwing 
    double amplitude_min = Amplitude_min;
    double amplitude_max = Amplitude_max;
    double sigma_min = Sigma_min;
    double sigma_max = Sigma_max;
    double scale_factor = 0;

    std::vector<GaussianKernel> Kernels;
    std::vector<int> surface_kernels;
    std::vector<Eigen::Vector3d> cached_pore_centers;
    std::vector<double> cached_pore_sdfs;
    std::vector<KernelRandomSpec> kernel_random_specs;
    std::vector<std::vector<int>> Paths;
    vector<pair<int, int>> max_paths_kernel;  //每个kernel通透性最大的路径两端
    vector<double> kernel_translucency;

    std::vector<Edge> Tube_edges;
    std::vector<std::vector<int>> Adj_list;
    std::vector<std::vector<int>> Unused_adj_list;

    //show compare
    std::vector<pair<int, int>> edge_con_mst;
    std::vector<pair<int, int>> edge_con_mbdst;
    std::vector<pair<int, int>> edge_con_msc;
    std::vector<pair<int, int>> edge_con_final;

    int model_solid_num = 0;
    int floating_voxel_removed = -1;
    double initPorosity = 0;
    double finalPorosity = 0;
	double finalTranslucency = 0;

    std::string outputPrefix = "D:/VSprojects/TaihuStone/result/" + input_file + "_" + std::to_string(PoresNum) + "_" + to_string_pre(Trans_thres, 2) + "_opt/";



};
