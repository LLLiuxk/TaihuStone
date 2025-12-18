#include "globalPara.h" 


int Resolution =100;

std::vector<double> Weights = { 1.0,1.0, 1.0 };  //mst 分类系数权重：边界-内部，内部-内部，边界-边界
std::vector<double> KT_weights = { 0.7, 0.1, 0.1, 0.1 };  //kernel translucency权重：角度、长度、内外分布、水平垂直
double Isolevel = 0;
double Gauss_level = 0.5;

int PoresNum = 30;
double surface_ratio = 0.3; //表面核占比

int Min_degree = 4;

double Amplitude_min = 1.0;
double Amplitude_max = 1.0;
double Sigma_min = 0.025;
double Sigma_max = 0.045;


double SmoothT = 10;    //控制平滑布尔的平滑区域大小, 越大，平滑效果越小，趋近于普通并集
double Tube_radius_factor = 0.4; // 控制管道半径相对于高斯核半径的比例
double Safe_distance_ratio = 0.3;

double Trans_thres = 0.8;  //kernel translucency threshold
double Adj_dis_thres = 0.3;
bool debug_show = false;
bool standard_show = true;


bool Direct_dis = false;

bool optimize_debug = true;