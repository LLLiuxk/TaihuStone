#include "globalPara.h" 

std::string input_file = "RockSetCr";
int Resolution = 100;

std::vector<double> Weights = { 0.8, 1.0, 1.2 };  //mst 分类系数权重：边界-内部，内部-内部，边界-边界
std::vector<double> KT_weights = { 0.7, 0.1, 0.1, 0.1 };  //kernel translucency权重：角度、长度、内外分布、水平垂直
double Isolevel = 0;
double Gauss_level = 0.5;

int PoresNum = 25;
double surface_ratio = 0.4; //表面核占比

int Min_degree = 4;

double Amplitude_min = 1.0;
double Amplitude_max = 1.0;
double Sigma_min = 0.015;
double Sigma_max = 0.033;
//250: 0.02, 0.04
//150~200: 0.025, 0.045
// 100: 0.033, 0.05
//120pores 0.02 0.037 0.58 0.38 A
//200pores 0.015 0.033 0.5 0.38 A
// 100 2 4 0.68 0.4 high_rock
//经验公式：sigma =(v/(6.8371*n*w))^(1/3)
// w for max = 5.0     for min = 30

double SmoothT = 15;    //控制平滑布尔的平滑区域大小, 越大，平滑效果越小，趋近于普通并集
double Tube_radius_factor = 0.9; // 控制管道半径相对于高斯核半径的比例, mid_radius_factor越大，通道越细
double Safe_distance_ratio = 0.7;

double Trans_thres = 0.8;  //kernel translucency threshold 
double Adj_dis_thres = 0.25;
bool debug_show = false;
bool standard_show = false;


bool Direct_dis = false;

bool optimize_debug = false;
bool Enable_noise = false;