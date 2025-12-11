#include "selfSupVoxel.h"

SupportCheckResult checkSupportVoxel(const VoxelGrid& grid, double densityThreshold = 0.5)
{
    SupportCheckResult result = { true, 0, 0.0, 0.0, 0.0 };

    // 1. 统计总体积
    long long totalSolidVoxels = 0;

    // 2. 遍历网格
    // 注意：k 从 0 开始遍历。k=0 的层我们也视为实体，但它贴在底板上，默认是被支撑的。
    // 所以检测逻辑主要针对 k >= 1 的层。
    for (int k = 0; k < grid.nz; ++k) {
        for (int j = 0; j < grid.ny; ++j) {
            for (int i = 0; i < grid.nx; ++i) {

                double val = grid.rho[grid.index(i, j, k)];

                // 如果当前格子是空的，跳过
                if (val < densityThreshold) continue;

                totalSolidVoxels++;

                // 如果是底层 (k=0)，默认由打印平台支撑，不需要额外支撑
                if (k == 0) continue;

                // 3. 检查下方支撑 (k-1 层)
                // 检查 3x3 邻域 (模拟 45 度悬垂角)
                // 如果 dx, dy, dz 差异很大，这里的逻辑可能需要调整范围
                bool isSupported = false;

                // 遍历下方 3x3 区域
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        // 获取下方体素的密度
                        // getSafe 处理了边界问题，所以不用担心 i+dx 越界
                        double belowVal = grid.getSafe(i + dx, j + dy, k - 1);

                        if (belowVal >= densityThreshold) {
                            isSupported = true;
                            goto check_done; // 只要找到一个支撑点，就跳出内层循环
                        }
                    }
                }

            check_done:;

                // 4. 判定结果
                if (!isSupported) {
                    result.unsupportedVoxelCount++;
                }
            }
        }
    }

    // 5. 计算统计数据
    double vVol = grid.voxelVolume();
    result.totalVolume = totalSolidVoxels * vVol;
    result.unsupportedVolume = result.unsupportedVoxelCount * vVol;

    if (result.totalVolume > 1e-9) {
        result.unsupportedRatio = result.unsupportedVolume / result.totalVolume;
    }
    else {
        result.unsupportedRatio = 0.0;
    }

    // 6. 最终判定
    // 设定一个容忍度，比如悬垂体积占比小于 0.1% 可能只是噪点或极小的倒角
    if (result.unsupportedRatio > 0.001) {
        result.isSupportFree = false;
    }

    return result;
}