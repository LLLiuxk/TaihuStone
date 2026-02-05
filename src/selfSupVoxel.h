#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <queue>
#include <Eigen/Dense>


struct VoxelGrid
{
    int nx, ny, nz;
    double dx, dy, dz;
    Eigen::Vector3d origin;
    std::vector<double> rho;

    VoxelGrid()
    {
        nx = ny = nz = 0;
        dx = dy = dz = 1.0;
		origin = Eigen::Vector3d::Zero();
    }
    VoxelGrid(int nx_, int ny_, int nz_,
        double dx_ = 1.0,
        double dy_ = 1.0,
        double dz_ = 1.0,
        const Eigen::Vector3d& origin_ = Eigen::Vector3d::Zero())
        : nx(nx_), ny(ny_), nz(nz_),
        dx(dx_), dy(dy_), dz(dz_),
        origin(origin_)
    {}
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

	//通过索引获取及修改体素密度
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
    bool isSupportFree = true;             // 是否免支撑
    int solid_nums = 0;
    int unsupportedVoxelCount = 0;// 需要支撑的体素数量
    double unsupportedRatio = 0.0;        // 悬垂比例
    int component_num = 0;
    std::vector<int> parts_solid_nums;
};

SupportCheckResult check_result_voxel(VoxelGrid& grid, double densityThreshold = 0.5);

SupportCheckResult checkSupportVoxel(VoxelGrid& grid, double densityThreshold = 0.5);

int findBaseLayer(const VoxelGrid& grid, double densityThreshold);

SupportCheckResult checkFloatingVoxel(VoxelGrid& grid, double densityThreshold = 0.5);



int removeFloatingSDF(Eigen::VectorXd& SDF, int nx, int ny, int nz, double smooth_radius, int& eliminate_num);

void getNeighbors26(int x, int y, int z, int nx, int ny, int nz, std::vector<Eigen::Vector3i>& nbrs);