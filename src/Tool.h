#pragma once
#include <igl/readOBJ.h>
#include <igl/signed_distance.h>

#include <igl/marching_cubes.h>
#include <igl/read_triangle_mesh.h>
#include <igl/voxel_grid.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/write_triangle_mesh.h>

#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <string>

#include "globalPara.h" 

using namespace std;
using namespace Eigen;

void Mesh2SDF(
    Eigen::MatrixXd& V,     // 网格顶点 (m×3)
    Eigen::MatrixXi& F,     // 网格面 (f×3)
    Eigen::MatrixXd& GV,    // 查询点集 (n×3)
    Eigen::VectorXd& S            // 输出: signed distance (n×1)
);


// SDF平滑并集。k越大，平滑效果越小，趋近于普通并集
double smoothUnionSDF(double sdf1, double sdf2, double k);

// SDF平滑交集
double smoothIntersecSDF(double sdf1, double sdf2, double k);

// SDF布尔运算：并集
double unionSDF(double sdf1, double sdf2);

// SDF布尔运算：交集
double intersectionSDF(double sdf1, double sdf2);

void view_model(Eigen::MatrixXd V1, Eigen::MatrixXi F1);
void view_two_models(Eigen::MatrixXd V1, Eigen::MatrixXi F1, Eigen::MatrixXd V2, Eigen::MatrixXi F2, Eigen::RowVector3d shift = RowVector3d(0.0, 0.0, 0.0));


void MarchingCubes(Eigen::VectorXd& S, Eigen::MatrixXd& GV, int nx, int ny, int nz, double isovalue, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

// SDF布尔运算：差集 (A - B)
double differenceSDF(double sdf1, double sdf2);

void align_models(const Eigen::MatrixXd& V1, const Eigen::MatrixXd& V2, Eigen::MatrixXd& V1_aligned, bool with_scale = false);

bool align_models_with_pca(const std::string& model1_path, const std::string& model2_path, const std::string& output_path);