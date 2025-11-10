#include "modelGen.h"

double GaussianKernel::gaussian_fun(const Eigen::Vector3d& p)
{
	double d2 = (p - center).squaredNorm();
	double gau_value = amplitude * std::exp(-d2 / (2.0 * sigma * sigma));   //sigma越大，函数越平缓
    //std::cout << "  gau_value:"<< gau_value <<std::endl;
	return gau_value;

}

GaussianKernel::GaussianKernel(Eigen::Vector3d center_, double sigma_, double amplitude_)
{
    center = center_;
    sigma = sigma_;         // 核的大小/影响力范围 (高斯函数的标准差)
    amplitude = amplitude_;
}

ModelGenerator::ModelGenerator(std::string input_file, int pores)
{
    if (!igl::read_triangle_mesh(input_file, V_ini, F_ini)) {
        std::cerr << "Error: Could not load model A." << std::endl;
        return;
    }
    pore_num = pores;
    //std::cout << "Model A loaded successfully." << std::endl;
    Mesh2SDF(V_ini, F_ini, GV, SDF_ini);
    generateGaussianSDF();
}


void ModelGenerator::generateGaussianSDF()
{
	// 随机生成空洞中心位置，pores = 0 采用不可预测随机数
	int pores = pore_num;
    pores = pores ? pores : (std::random_device{}()%PoresNum);
    std::cout << "pores: " << pores << std::endl;
    std::mt19937 gen(pores);

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

	double model_count = inside_indices.size();

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
        Eigen::Vector3d candidate_center = GV.row(chosen_idx).transpose();

        // 检查与已有空洞中心的最小距离
        bool valid = true;
        Eigen::Vector3d min_pt = GV.colwise().minCoeff();
        Eigen::Vector3d max_pt = GV.colwise().maxCoeff();
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

    std::cout << "Generate " << pore_centers.size() << " kernels" << std::endl;
    //for (int i = 0; i < pore_centers.size(); i++)   std::cout << "i: " << i << "  " << pore_centers[i] << std::endl;

    std::uniform_real_distribution<double> dist_amp(amplitude_min, amplitude_max);
    std::uniform_real_distribution<double> dist_sigma(sigma_min, sigma_max);

    // 为每个空洞中心生成随机参数
    std::vector<double> pore_amplitudes;
    std::vector<double> pore_sigmas;

    int pore_size = pore_centers.size();
    for (size_t i = 0; i < pore_size; ++i) {
		double sigma_val = dist_sigma(gen);
		double amplitude_val = dist_amp(gen);
        std::cout << "i: " << i << "  " << sigma_val<<"  "<<amplitude_val << std::endl;
        GaussianKernel kernel(pore_centers[i], sigma_val, amplitude_val);
        Kernels.push_back(kernel);
        /*pore_amplitudes.push_back(dist_amp(gen));
        pore_sigmas.push_back(dist_sigma(gen));*/
    }

    std::cout << "Combine the Gaussian fileds... " << std::endl;
    double void_count = 0;
    SDF_out.resize(grid_num);
	//单独保存高斯孔隙场SDF
    Eigen::VectorXd SDF_gaussian(grid_num);
    for (int idx = 0; idx < grid_num; ++idx) {
        Eigen::Vector3d p = GV.row(idx);
		SDF_gaussian(idx) = combinedSDF(p);
        SDF_out(idx) = smooth_IntersecSDF(SDF_ini(idx), -SDF_gaussian(idx), smooth_t);
        //SDF_out(idx) = differenceSDF(SDF_ini(idx), SDF_gaussian(idx));
        //SDF_out(idx) = differenceSDF(SDF_ini(idx), 0);
        if (SDF_gaussian(idx) < isolevel)
           // std::cout << "idx:  " << idx << "   SDF_ini(idx):" << SDF_ini(idx) << "   SDF_gaussian(idx):  " << SDF_gaussian(idx) <<"  SDF_out(idx) :" << SDF_out(idx) << std::endl;
        //std::cout << SDF_out(idx) << std::endl;
        if (SDF_out(idx) < isolevel) {
            void_count += 1;
        }
    }
    //std::cout << "成功在仿生形状内放置 " << void_centers.size() << " 个空洞点" << std::endl;
 
     // 计算孔隙率
    double porosity = void_count / model_count;
    std::cout << "Porosity: " << porosity * 100 << "%" << std::endl;

    // 存储最终的孔隙率
    finalPorosity = porosity;

    // Marching Cubes
    MarchingCubes(SDF_out, GV, resolution, resolution, resolution, isolevel, V_out, F_out);

    Eigen::MatrixXd V_g; //输出网格顶点
    Eigen::MatrixXi F_g; // 输出网格面片
    MarchingCubes(SDF_gaussian, GV, resolution, resolution, resolution, isolevel, V_g, F_g);
    int comp = single_component(V_g, F_g);
    //view_model(V_g, F_g);
	//std::string filename = "output/gaussian_pores.stl";
 //   std::filesystem::path filePath(filename);
 //   std::filesystem::path dir = filePath.parent_path();
 //   if (!dir.empty() && !std::filesystem::exists(dir)) {
 //       std::filesystem::create_directories(dir);
 //   }

    // 保存STL文件
    //bool stl_success = igl::writeSTL(filename, V, F, N, igl::FileEncoding::Binary);

    view_two_models(V_out, F_out, V_g, F_g, Eigen::RowVector3d(1, 0,0));
    // calculate normal
    //Eigen::MatrixXd N;
    //igl::per_face_normals(V_out, F_out, N);
    std::cout << "Generated mesh: " << V_out.rows() << " vertices, " << F_out.rows() << " faces" << std::endl;
    return;


}

double ModelGenerator::combinedSDF(Eigen::Vector3d & p)
{
    double total_void = 0.0;
    int gaus_num = Kernels.size();
    for (size_t i = 0; i < gaus_num; i++) {
        //std::cout <<"gaus_num: "<< gaus_num<<"  "<< i << std::endl;
        total_void += Kernels[i].gaussian_fun(p);
    }
    //std::cout << "total_void: " << total_void << std::endl;
    return  gauss_combined - total_void;  // 当前使用：负的空洞总和 
}

void ModelGenerator::show_model()
{
    view_two_models(V_ini, F_ini, V_out, F_out, Eigen::RowVector3d(1, 0.0, 0.0));

}