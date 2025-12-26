#include "selfSupVoxel.h"

SupportCheckResult checkSupportVoxel(VoxelGrid& grid, double densityThreshold)
{
    SupportCheckResult result = { true, 0, 0.0};

    // 1. 统计总体积
    int totalSolidVoxels = 0;

    // 2. 遍历网格
    // 注意：k 从 0 开始遍历。k=0 的层我们也视为实体，但它贴在底板上，默认是被支撑的。
    std::vector<bool> supported(grid.nx * grid.ny * grid.nz, false);
	int base_layer = findBaseLayer(grid, densityThreshold);
   // for (int j = 0; j < grid.ny; ++j) {
   //     for (int i = 0; i < grid.nx; ++i) {
			//int idx = grid.index(i, j, base_layer);
   //         supported[idx] = true;
   //     }
   // }

    for (int k = base_layer; k < grid.nz; ++k) {
        for (int j = 0; j < grid.ny; ++j) {
            for (int i = 0; i < grid.nx; ++i) {
                int idx = grid.index(i, j, k);
                double val = grid.rho[grid.index(i, j, k)];

                // 如果当前格子是空的，跳过
                if (val < densityThreshold) continue;

                totalSolidVoxels++;

                // 如果是底层 (k=0)，默认由打印平台支撑，不需要额外支撑
                if (k == base_layer) 
                {
                    supported[idx] = true;
                    continue;
                }

                // 3. 检查下方支撑 (k-1 层)， 检查 3x3 邻域 (模拟 45 度悬垂角)

                // 遍历下方 3x3 区域
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        // 获取下方体素的密度
                        double belowVal = grid.getSafe(i + dx, j + dy, k - 1);
                        if (belowVal >= densityThreshold) { // && supported[grid.index(i + dx, j + dy, k - 1)]
                            supported[idx] = true;
                            goto check_done; // 只要找到一个支撑点，就跳出内层循环
                        }
                    }
                }
            check_done:;
                if (!supported[idx]) {
                    result.unsupportedVoxelCount++;
                    //grid.at(i, j, k) = 0.0;
                }
            }
        }
    }
    int usnum = 0;
    for(auto p: supported)
        if(p)
            usnum++;

	std::cout << "unsupport num: " << totalSolidVoxels - usnum << "   " << result.unsupportedVoxelCount << std::endl;
    // 5. 计算统计数据
    //double vVol = grid.voxelVolume();
    //result.totalVolume = totalSolidVoxels * vVol;
    //result.unsupportedVolume = result.unsupportedVoxelCount * vVol;
    //if (result.totalVolume > 1e-9) {
    //    result.unsupportedRatio = result.unsupportedVolume / result.totalVolume;
    //}
    //else {
    //    result.unsupportedRatio = 0.0;
    //}
	//std::cout << "result.unsupportedVoxelCount" << result.unsupportedVoxelCount << "   totalSolidVoxels: " << totalSolidVoxels << std::endl;
    result.unsupportedRatio = 1.0 * result.unsupportedVoxelCount / totalSolidVoxels;
    // 6. 最终判定
    // 设定一个容忍度，比如悬垂体积占比小于 0.1% 可能只是噪点或极小的倒角
    if (result.unsupportedRatio > 0.001) {
        result.isSupportFree = false;
    }
    if (result.isSupportFree) std::cout << "Result is SupportFree!" << std::endl;
    else std::cout << "Result NEED Support!" << std::endl
        << result.unsupportedVoxelCount << "   " << result.unsupportedRatio *100<<"%"<< std::endl;

    return result;
}


int findBaseLayer(const VoxelGrid& grid, double densityThreshold)
{
    for (int k = 0; k < grid.nz; ++k) {
        for (int j = 0; j < grid.ny; ++j) {
            for (int i = 0; i < grid.nx; ++i) {
                double val = grid.rho[grid.index(i, j, k)];
                // 如果当前格子是空的，跳过
                if (val < densityThreshold) continue;
                // 找到第一层实体体素
                else
                    return k;
            }
        }
    }
	return -1; // 未找到实体体素
}