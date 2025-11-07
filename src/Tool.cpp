#include "Tool.h"


void Mesh2SDF(Eigen::MatrixXd& V,  Eigen::MatrixXi& F, Eigen::MatrixXd& GV, Eigen::VectorXd& SDF)
{
    if (V.rows() == 0 || F.rows() == 0 ) {
        std::cerr << " Input failed! Wrong mesh!" << std::endl;
        return;
    }

    Eigen::VectorXi I;  // 最近面索引
    Eigen::MatrixXd C;  // 最近点坐标
    Eigen::MatrixXd N;  // 内外法向符号

    Eigen::Vector3d bb_min = V.colwise().minCoeff();
    Eigen::Vector3d bb_max = V.colwise().maxCoeff();
    Eigen::Vector3d bb_size = bb_max - bb_min;
    
    //normalize: scaling and moving
    double max_dim = bb_size.maxCoeff();
    double scale_factor = 1.0 / (max_dim + 1e-9);
    Eigen::Vector3d bb_center = (bb_min + bb_max) / 2.0;
    V = (V.rowwise() - bb_center.transpose()) * scale_factor;
    bb_min = V.colwise().minCoeff();

    int res = Resolution;  // 网格分辨率，可调高以提升精度
    double dx = 1.0 / (res - 1);

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

     std::cout << "Sample resolution: " << res<<"  "<< res<<"  "<< res<<" = "<< GV.rows() << std::endl;

    // 调用 libigl 的 signed_distance()
    igl::signed_distance( GV, V, F, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, SDF, I, C, N );
}


// SDF平滑并集。k越大，平滑效果越小，趋近于普通并集
double smoothUnionSDF(double sdf1, double sdf2, double k)
{
    return -log(exp(-sdf1 * k) + exp(-sdf2 * k)) / k;
}

// SDF平滑交集
double smoothIntersecSDF(double sdf1, double sdf2, double k)
{
    return  log(exp(-sdf1 * k) + exp(-sdf2 * k)) / k;
}
