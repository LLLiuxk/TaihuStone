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



