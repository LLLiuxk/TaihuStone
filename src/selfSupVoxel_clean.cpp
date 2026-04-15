#include "selfSupVoxel.h"

SupportCheckResult checkSupportVoxel(VoxelGrid& grid, double densityThreshold)
{
    SupportCheckResult result;

    // 1.  
    int totalSolidVoxels = 0;

    // 2. 
    //  k  0  k=0  
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

                //  
                if (val < densityThreshold) continue;

                totalSolidVoxels++;

                //   (k=0) 
                if (k == base_layer) 
                {
                    supported[idx] = true;
                    continue;
                }

                // 3.   (k-1 )  3x3  (  45 )

                //   3x3 
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        //  
                        double belowVal = grid.getSafe(i + dx, j + dy, k - 1);
                        if (belowVal >= densityThreshold /*&& supported[grid.index(i + dx, j + dy, k - 1)]*/) { // && supported[grid.index(i + dx, j + dy, k - 1)]
                            supported[idx] = true;
                            goto check_done; //  
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

	result.solid_nums = totalSolidVoxels;
	std::cout << "unsupport num: " << totalSolidVoxels<<" - " <<usnum << " =  " << result.unsupportedVoxelCount << std::endl;
    // 5.  
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
    // 6.  
    //   0.1%  
    if (result.unsupportedRatio > 0.001) {
        result.isSupportFree = false;
    }
    if (result.isSupportFree) std::cout << "Result is SupportFree!" << std::endl;
    else std::cout << "Result NEED Support!" << std::endl;

    std::cout  << "Unsupport voxel num: "<<result.unsupportedVoxelCount << "   ratio: " << result.unsupportedRatio * 100 << "%" << std::endl;

    return result;
}


int findBaseLayer(const VoxelGrid& grid, double densityThreshold)
{
    for (int k = 0; k < grid.nz; ++k) {
        for (int j = 0; j < grid.ny; ++j) {
            for (int i = 0; i < grid.nx; ++i) {
                double val = grid.rho[grid.index(i, j, k)];
                //  
                if (val < densityThreshold) continue;
                //  
                else
                    return k;
            }
        }
    }
	return -1; //  
}

SupportCheckResult checkFloatingVoxel(VoxelGrid& grid, double densityThreshold)
{
    SupportCheckResult result;
    const int nx = grid.nx;
    const int ny = grid.ny;
    const int nz = grid.nz;

    const int N = nx * ny * nz;

    std::vector<int> parts_solids;
    

    // visited 
    std::vector<char> visited(N, 0);
    // 6 
    const int dx[6] = { 1,-1, 0, 0, 0, 0 };
    const int dy[6] = { 0, 0, 1,-1, 0, 0 };
    const int dz[6] = { 0, 0, 0, 0, 1,-1 };

    auto inside = [&](int i, int j, int k) {
        return (i >= 0 && i < nx &&
                j >= 0 && j < ny &&
                k >= 0 && k < nz);
    };

    // Flood fill
    for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
        int id = grid.index(i, j, k);

        //  
        if (visited[id]) continue;
        if (grid.at(i, j, k) <= densityThreshold) continue;

        int component_count = 0;
        std::queue<Eigen::Vector3i> q;

        q.emplace(i, j, k);
        visited[id] = 1;

        while (!q.empty())
        {
            Eigen::Vector3i p = q.front();
            q.pop();
            component_count++;

            for (int d = 0; d < 6; ++d)
            {
                int ni = p.x() + dx[d];
                int nj = p.y() + dy[d];
                int nk = p.z() + dz[d];

                if (!inside(ni, nj, nk)) continue;

                int nid = grid.index(ni, nj, nk);
                if (visited[nid]) continue;
                if (grid.at(ni, nj, nk) <= densityThreshold) continue;

                visited[nid] = 1;
                q.emplace(ni, nj, nk);
            }
        }

        parts_solids.push_back(component_count);
    }

    //  
    std::sort(parts_solids.begin(),
        parts_solids.end(),
        std::greater<int>());

    result.component_num = parts_solids.size();
	result.parts_solid_nums = parts_solids;
    int floating_nums = 0;
    for (auto ps : parts_solids)
        floating_nums += ps;
    std::cout << "This result has " << result.component_num << " components, and the biggest one has " << parts_solids[0] << " voxels, other floating voxels are: " << floating_nums - parts_solids[0] << " voxels"<<std::endl;
    return result;
}


SupportCheckResult check_result_voxel(VoxelGrid& grid, double densityThreshold)
{
    SupportCheckResult result = checkSupportVoxel(grid, densityThreshold);

    SupportCheckResult result2 = checkFloatingVoxel(grid, densityThreshold);

	result.component_num = result2.component_num;
	result.parts_solid_nums = result2.parts_solid_nums;
     
    return result;
}


int removeFloatingSDF(Eigen::VectorXd& SDF, int nx, int ny, int nz, double smooth_radius, int& eliminate_num)
{
    int N = nx * ny * nz;
    std::vector<int> comp_id(N, -1);
    std::vector<int> comp_size;
    int cid = 0;

    std::vector<Eigen::Vector3i> stack, nbrs;

    // ===== 1.  SDF <= 0=====
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
            {
                int id = x + nx * (y + ny * z); //idx(x, y, z, g);
                if (SDF[id] <= 0 && comp_id[id] == -1)
                {
                    int count = 0;
                    stack.clear();
                    stack.emplace_back(x, y, z);
                    comp_id[id] = cid;

                    while (!stack.empty())
                    {
                        auto p = stack.back();
                        stack.pop_back();
                        ++count;

                        getNeighbors26(p.x(), p.y(), p.z(), nx, ny, nz, nbrs);
                        for (auto& q : nbrs)
                        {
                            int qid = q.x() + nx * (q.y() + ny * q.z());
                            if (SDF[qid] <= 0 && comp_id[qid] == -1)
                            {
                                comp_id[qid] = cid;
                                stack.push_back(q);
                            }
                        }
                    }
                    comp_size.push_back(count);
                    cid++;
                }
            }

    // ===== 2.   =====
    int main_comp = 0;
    for (int i = 1; i < cid; ++i)
        if (comp_size[i] > comp_size[main_comp])
            main_comp = i;

    // ===== 3.   =====
    std::vector<int> main_voxels;
    for (int i = 0; i < N; ++i)
        if (comp_id[i] == main_comp)
            main_voxels.push_back(i);

    // ===== 4.   =====
    int modified = 0;
    int removed_components = cid - 1;

    for (int i = 0; i < N; ++i)
    {
        if (comp_id[i] >= 0 && comp_id[i] != main_comp)
        {
            //   O(N)  KD-tree
            double min_d2 = 1e30;

            int z = i / (nx * ny);
            int y = (i - z * nx * ny) / nx;
            int x = i % nx;

            for (int j : main_voxels)
            {
                int mz = j / (nx * ny);
                int my = (j - mz * nx * ny) / nx;
                int mx = j % nx;

                double dx = (x - mx);
                double dy = (y - my);
                double dz = (z - mz);
                min_d2 = std::min(min_d2, dx * dx + dy * dy + dz * dz);
            }

            double d = std::sqrt(min_d2);
            double t = std::min(d / smooth_radius, 1.0);

            // smoothstep  
            double s = t * t * (3 - 2 * t);
            SDF[i] = s * std::abs(SDF[i]);

            modified++;
        }
    }
    eliminate_num = modified;
    return removed_components;
}

void getNeighbors26(int x, int y, int z, int nx, int ny, int nz, std::vector<Eigen::Vector3i>& nbrs)
{
    nbrs.clear();
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz)
            {
                if (dx == 0 && dy == 0 && dz == 0) continue;
                int nx_ = x + dx, ny_ = y + dy, nz_ = z + dz;
                if (nx_ >= 0 && nx_ < nx &&
                    ny_ >= 0 && ny_ < ny &&
                    nz_ >= 0 && nz_ < nz)
                {
                    nbrs.emplace_back(nx_, ny_, nz_);
                }
            }
}