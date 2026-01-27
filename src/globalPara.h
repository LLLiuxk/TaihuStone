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
extern double surface_ratio;

extern int Min_degree;
extern double Amplitude_min;
extern double Amplitude_max;
extern double Sigma_min;
extern double Sigma_max;
extern double W4sig_max; //7.0
extern double W4sig_min; // 42.0;
extern double Axis_max_ratio;


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
extern bool Direct_dis;
extern bool Iso_kernel;  //使用iso还是椭球
extern bool Handle_overlap;
extern bool optimize_debug;
extern bool Enable_noise;
extern bool topo_optimize;
extern bool dynamic_change_para;


std::string trim(const std::string& str);

// 辅助函数：解析 vector<double> (例如 "0.8, 1.0, 1.2")
std::vector<double> parseVector(std::string valStr);

bool loadParameters(const std::string& filename);
