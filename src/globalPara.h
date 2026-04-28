#pragma once// 参数结构体
#define NO_DEBUG false
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

extern std::string input_file;

extern int Resolution;

extern std::vector<double> Weights;
extern std::vector<double> KT_weights;
extern double Isolevel;

extern int PoresNum;
extern int Max_edge_limit;
extern double surface_ratio;

extern int Max_degree;
extern double Amplitude_min;
extern double Amplitude_max;
extern double Sigma_min;
extern double Sigma_max;
extern double W4sig_max; //7.0
extern double W4sig_min; // 42.0;
extern double Axis_max_ratio;
extern double BaseLayer;

extern double Gauss_level;

extern double SmoothT;    //控制平滑布尔的平滑区域大小

extern double Tube_radius_factor;  // 控制管道半径相对于高斯核半径的比例

extern double Safe_distance_ratio; // Dart throwing安全距离比例

extern double Trans_thres;

extern double Adj_dis_thres;
extern bool debug_show;
extern bool standard_show;
extern bool figure_show;
extern bool compare_show;
extern bool scr_figure;
extern bool Direct_dis;
extern bool Iso_kernel;  //使用iso还是椭球
extern bool Handle_overlap;
extern bool optimize_debug;
extern bool Enable_noise;
extern bool topo_optimize;
extern bool dynamic_change_para;
extern bool fixed_graph_dual_render;
extern bool auto_capture_translucency_graph;

extern double Param2_Amplitude_min;
extern double Param2_Amplitude_max;
extern double Param2_Sigma_min;
extern double Param2_Sigma_max;
extern double Param2_W4sig_max;
extern double Param2_W4sig_min;
extern double Param2_Axis_max_ratio;
extern double Param2_Gauss_level;
extern double Param2_SmoothT;
extern double Param2_Tube_radius_factor;

extern float show_degree_x;
extern float show_degree_y;
extern double show_degree_z;

std::string trim(const std::string& str);

// 辅助函数：解析 vector<double> (例如 "0.8, 1.0, 1.2")
std::vector<double> parseVector(std::string valStr);

bool loadParameters(const std::string& filename);

void writeVector(std::ofstream& ofs, const std::string& name, const std::vector<double>& vec);

void saveParameters(const std::string& filename);

inline int b2i(bool v) { return v ? 1 : 0; }

std::string generateFilename();
