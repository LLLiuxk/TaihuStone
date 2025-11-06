#pragma once
#include <igl/readOBJ.h>
#include <igl/signed_distance.h>

#include <igl/marching_cubes.h>
#include <igl/read_triangle_mesh.h>
#include <igl/voxel_grid.h>
#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <string>

#include "globalPara.h" 


void Mesh2SDF(
    const Eigen::MatrixXd& V,     // 网格顶点 (m×3)
    const Eigen::MatrixXi& F,     // 网格面 (f×3)
    Eigen::MatrixXd& GV,    // 查询点集 (n×3)
    Eigen::VectorXd& S            // 输出: signed distance (n×1)
);

// SDF平滑并集。k越大，平滑效果越小，趋近于普通并集
double smoothUnionSDF(double sdf1, double sdf2, double k);

// SDF平滑交集
double smoothIntersecSDF(double sdf1, double sdf2, double k);

double ModelGenerator::combinedSDF(const Eigen::RowVector3d& p,
    const Eigen::Vector3d& radii,
    const std::vector<Eigen::Vector3d>& void_centers,
    const std::vector<double>& void_amplitudes,
    const std::vector<double>& void_sigmas,
    const double x,
    const double t) const {
    double total_void = 0.0;
    for (size_t i = 0; i < void_centers.size(); i++) {
        total_void += gaussianKernel(p, void_centers[i], void_amplitudes[i], void_sigmas[i]);
    }

