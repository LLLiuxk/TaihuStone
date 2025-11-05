#include "Tool.h"


void Mesh2SDF(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXd& GV, Eigen::VectorXd& SDF)
{
    if (V.rows() == 0 || F.rows() == 0 ) {
        std::cerr << "[OBJ2SDF] 输入网格或采样点为空！" << std::endl;
        return;
    }

    Eigen::VectorXi I;  // 最近面索引
    Eigen::MatrixXd C;  // 最近点坐标
    Eigen::MatrixXd N;  // 内外法向符号

    Eigen::Vector3d bb_min = V.colwise().minCoeff();
    Eigen::Vector3d bb_max = V.colwise().maxCoeff();
    Eigen::Vector3d bb_size = bb_max - bb_min;

    int res = Resolution;  // 网格分辨率，可调高以提升精度
    double dx = bb_size.maxCoeff() / (res - 1);

	int total_points = res * res * res;
    GV.resize(total_points, 3);
    SDF.resize(total_points);
    /*Eigen::MatrixXd GV(total_points, 3);
    Eigen::VectorXd SDF(total_points);*/
    double void_count = 0;
 
   std:: vector<Eigen::Vector3d> grid_points;
    grid_points.reserve(res * res * res);
        for (int i = 0; i < res; ++i) {
        for (int j = 0; j < res; ++j) {
            for (int k = 0; k < res; ++k) {
                Eigen::Vector3d p = bb_min + Eigen::Vector3d(i, j, k) * dx;
                grid_points.push_back(p);
            }
        }
    }
    for (int i = 0; i < grid_points.size(); ++i)
        GV.row(i) = grid_points[i];

     std::cout << "采样点数量: " <<  GV.rows() << std::endl;

    // 调用 libigl 的 signed_distance()
    igl::signed_distance( GV, V, F, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, SDF, I, C, N );
}