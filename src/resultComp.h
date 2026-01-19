#pragma once
#ifndef RESULT_COMP_H
#define RESULT_COMP_H

#include <vector>
#include <set>
#include <Eigen/Core>

// 简单的 3D 向量结构，用于传入物理坐标范围
struct Vec3 {
    double x, y, z;
    double distSq(const Vec3& other) const {
        return (x - other.x) * (x - other.x) + (y - other.y) * (y - other.y) + (z - other.z) * (z - other.z);
    }
};

class GaussianFieldMSC {
public:
    /**
     * @brief 构造函数
     * @param nx X轴网格分辨率
     * @param ny Y轴网格分辨率
     * @param nz Z轴网格分辨率
     * @param min_b 物理空间的最小包围盒坐标
     * @param max_b 物理空间的最大包围盒坐标
     */
    GaussianFieldMSC(int nx, int ny, int nz, Vec3 min_b, Vec3 max_b);

    /**
     * @brief 核心函数：计算 MSC 连接关系
     *
     * @param field_data 展平的标量场数据 (大小应为 nx*ny*nz，顺序 X->Y->Z)
     * @param original_points 原始的高斯核中心点 (用于最后将拓扑结构映射回这些点)
     * @return std::vector<std::vector<int>> 邻接表，索引对应 original_points 的下标
     */
    std::vector<std::vector<int>> ComputeConnectivity(
        const Eigen::VectorXd& field_data,
        const std::vector<Vec3>& original_points);

private:
    int nx_, ny_, nz_;
    Vec3 min_b_, max_b_;
    double dx_, dy_, dz_;

    // 内部辅助函数
    inline int getIdx(int x, int y, int z) const;
    inline void getCoord(int idx, int& x, int& y, int& z) const;

    bool CheckIsMaxima(const Eigen::VectorXd& field, int x, int y, int z);
    bool CheckIs2Saddle(const Eigen::VectorXd& field, int x, int y, int z);

    std::set<int> TraceAscent(const Eigen::VectorXd& field, const std::vector<bool>& is_maxima, int start_idx);
    int MapGridToPoint(int grid_idx, const std::vector<Vec3>& points);
};

#endif // RESULT_COMP_H