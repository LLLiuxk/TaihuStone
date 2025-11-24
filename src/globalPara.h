#pragma on// 参数结构体

extern int Resolution;

extern double Isolevel;

extern int PoresNum;

extern int Min_degree;
extern double Amplitude_min;
extern double Amplitude_max;
extern double Sigma_min;
extern double Sigma_max;
extern double Gauss_level;

extern double SmoothT;    //控制平滑布尔的平滑区域大小

extern double Tube_radius_factor;  // 控制管道半径相对于高斯核半径的比例

extern double Safe_distance_ratio; // Dart throwing安全距离比例

extern bool debug_show;
extern bool Direct_dis;
//struct Parameters {
//    Eigen::Vector3d radii = { 1, 0.5, 0.5 };
//    int num_voids = 188;
//    double min_distance = 0.15;
//    double amplitude_min = 1;
//    double amplitude_max = 1;
//    double sigma_min = 150;
//    double sigma_max = 150;
//    int resolution = 100;
//    double bound = 1.2;
//    double isolevel = 0;
//    double t = 0.1;    //控制平滑布尔的平滑区域大小
//    unsigned int seed = 0;
//    double porosity = 0.7; // 孔隙率
//    double yta = 0.01;    // 误差
//    bool enable_self_support_post = false; // 是否启用自支撑后处理（默认关闭）
//
//    // 新增：基础形状参数
//    enum BaseShapeType {
//        ELLIPSOID,
//        CYLINDER,
//        SPHERE
//    };
//    BaseShapeType base_shape = ELLIPSOID;
//
//    // 表面采样参数
//    int surface_sample_count = 100;      // 表面采样点数量
//    double surface_gaussian_amplitude = 1;  // 表面高斯核振幅
//    double surface_gaussian_sigma = 100;     // 表面高斯核标准差
//}; 
