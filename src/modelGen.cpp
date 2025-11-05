#include "modelGen.h"

double GaussianKernel::gaussian_fun(const Eigen::RowVector3d& p)
{
	double d2 = (p - center.transpose()).squaredNorm();
	return amplitude * std::exp(-d2 * sigma);

}


void ModelGenerator::generateGaussianSDF(int pores)
{
	// 随机生成空洞中心位置，pores = 0 采用不可预测随机数
    pores = pores ? pores : (std::random_device{}()%PoresNum);
    std::cout << "pores: " << pores << std::endl;
    std::mt19937 gen(pores);
    std::uniform_real_distribution<double> dist_amp(amplitude_min, amplitude_max);
    std::uniform_real_distribution<double> dist_sigma(sigma_min, sigma_max);

    int grid_num = SDF_ini.size();
    if (grid_num == 0) {
        std::cerr << "Empty Initial SDF!" << std::endl;
        return;
	}
    // search inside points
    std::vector<Eigen::Vector3d> pore_centers;
    std::vector<int> inside_indices;
    for (int idx = 0; idx < grid_num; ++idx) {
        if (SDF_ini(idx) < isolevel) {
            inside_indices.push_back(idx);
        }
    }

    if (inside_indices.empty()) {
        std::cerr << "Warning: no legal points inside!" << std::endl;
        return;
    }

    // 在整个形状内部采样
    std::uniform_int_distribution<int> index_dist(0, inside_indices.size() - 1);

    int attempts = 0;
    const int max_attempts = pores * 50;

    while (pore_centers.size() < pores && attempts < max_attempts) {
        attempts++;

        // 随机选择一个内部点
        int chosen_idx = inside_indices[index_dist(gen)];
        Eigen::Vector3d candidate_center = GV_ini.row(chosen_idx).transpose();

        // 检查与已有空洞中心的最小距离
        bool valid = true;
        Eigen::Vector3d min_pt = GV_ini.colwise().minCoeff();
        Eigen::Vector3d max_pt = GV_ini.colwise().maxCoeff();
        Eigen::Vector3d box_size = max_pt - min_pt;
        double volume = box_size.x() * box_size.y() * box_size.z();
        safe_distance = 0.5 * std::cbrt(volume / pores);

        for (const auto& existing_center : pore_centers) {
            if ((candidate_center - existing_center).squaredNorm() <  safe_distance* safe_distance) {
                valid = false;
                break;
            }
        }

        if (valid) {
            pore_centers.push_back(candidate_center);
        }
    }

    std::cout << "Generate " << pore_centers.size() << " pores" << std::endl;
 
}