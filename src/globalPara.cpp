#include "globalPara.h" 


int Resolution = 100;

double Isolevel = 0.0;

int PoresNum = 100;

double Amplitude_min = 1.0;
double Amplitude_max = 1.2;
double Sigma_min = 0.05;
double Sigma_max = 0.1;
double Gauss_level = Amplitude_min * 0.5;

double SmoothT = 2;    //控制平滑布尔的平滑区域大小