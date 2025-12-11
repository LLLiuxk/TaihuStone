#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "Tool.h"

struct VoxelGrid
{
    int nx, ny, nz;
    double dx, dy, dz;
    Eigen::Vector3d origin;
    std::vector<double> rho;

    inline int index(int i, int j, int k) const
    {
        return i + nx * (j + ny * k);
    }

    // 为了安全起见，添加一个const版本的访问器，并带边界检查
    double getSafe(int i, int j, int k) const
    {
        if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz)
            return 0.0; // 越界视为空气
        return rho[index(i, j, k)];
    }

    double& at(int i, int j, int k)
    {
        return rho[index(i, j, k)];
    }

    // 获取单个体素的体积
    double voxelVolume() const {
        return dx * dy * dz;
    }
};

// 检测结果结构体
struct SupportCheckResult {
    bool isSupportFree;             // 是否免支撑
    long long unsupportedVoxelCount;// 需要支撑的体素数量
    double unsupportedVolume;       // 需要支撑的体积
    double totalVolume;             // 模型总体积
    double unsupportedRatio;        // 悬垂比例
};

SupportCheckResult checkSupportVoxel(const VoxelGrid& grid, double densityThreshold = 0.5); 