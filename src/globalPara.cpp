#include "globalPara.h" 


int Resolution =100;

double Isolevel = 0;
double Gauss_level = 0.5;

int PoresNum = 20;

int Min_degree = 4;

double Amplitude_min = 1.0;
double Amplitude_max = 1.0;
double Sigma_min = 0.025;
double Sigma_max = 0.04;


double SmoothT = 50;    //控制平滑布尔的平滑区域大小, 越大，平滑效果越小，趋近于普通并集
double Tube_radius_factor = 0.4; // 控制管道半径相对于高斯核半径的比例
double Safe_distance_ratio = 0.35;

double Trans_thres = 0.6;

bool debug_show = false;
bool standard_show = true;
bool Direct_dis = false;