#pragma once
#include <igl/readOBJ.h>
#include <igl/signed_distance.h>
#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <string>

#include "globalPara.h" 


void OBJ2SDF(
    const Eigen::MatrixXd& V,     // 网格顶点 (m×3)
    const Eigen::MatrixXi& F,     // 网格面 (f×3)
    Eigen::MatrixXd& GV,    // 查询点集 (n×3)
    Eigen::VectorXd& S            // 输出: signed distance (n×1)
);

