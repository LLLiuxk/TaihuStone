#pragma once
#include <Eigen/Dense>
#include <vector>


// 类型别名
using Vector3d = Eigen::Vector3d;

// 辅助函数声明
namespace MorseComplex {
    // 主要功能函数
    std::vector<std::pair<int, int>>  compare_msc(std::vector<Vector3d> Kernels,
        Eigen::VectorXd SDF_gaussian,
        int res,
        int grid_num);
}