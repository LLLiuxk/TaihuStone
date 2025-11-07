#include "modelGen.h"

double GaussianKernel::gaussian_fun(const Eigen::Vector3d& p)
{
	double d2 = (p - center).squaredNorm();
    std::cout << p << "  center:  " << center << "  d2:"<<d2<<std::endl;
	return amplitude * std::exp(-d2 /(2.0* sigma* sigma));

}

GaussianKernel::GaussianKernel(Eigen::Vector3d center_, double sigma_, double amplitude_)
{
    center = center_;
    sigma = sigma_;         // 核的大小/影响力范围 (高斯函数的标准差)
    amplitude = amplitude_;
}

ModelGenerator::ModelGenerator(std::string input_file)
{
    if (!igl::read_triangle_mesh(input_file, V_ini, F_ini)) {
        std::cerr << "Error: Could not load model A." << std::endl;
        return;
    }
    //std::cout << "Model A loaded successfully." << std::endl;
    Mesh2SDF(V_ini, F_ini, GV, SDF_ini);
    generateGaussianSDF();
}


void ModelGenerator::generateGaussianSDF(int pores)
{
    auto start_time = std::chrono::high_resolution_clock::now();

	// 随机生成空洞中心位置，pores = 0 采用不可预测随机数
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

	double ellipsoid_count = inside_indices.size();

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
		//GaussianKernel kernel(pore_centers[i], dist_sigma(gen), dist_amp(gen));
        Kernels.push_back(kernel);
        /*pore_amplitudes.push_back(dist_amp(gen));
        pore_sigmas.push_back(dist_sigma(gen));*/
    }

    std::cout << "Combine the Gaussian fileds... " << std::endl;
    double void_count = 0;
    for (int idx = 0; idx < grid_num; ++idx) {
        Eigen::Vector3d p = GV.row(idx);
        std::cout << "idx:  " << idx << std::endl;
        SDF_out(idx) = smoothIntersecSDF(SDF_ini(idx), -combinedSDF(p), smooth_t);
        if (SDF_out(idx) < isolevel) {
            void_count += 1;
        }
    }
    //std::cout << "成功在仿生形状内放置 " << void_centers.size() << " 个空洞点" << std::endl;
 
     // 计算孔隙率
    double porosity = void_count / ellipsoid_count;
    std::cout << "Porosity: " << porosity * 100 << "%" << std::endl;

    // 存储最终的孔隙率
    finalPorosity = porosity;

    // Marching Cubes
    igl::marching_cubes(SDF_out, GV, resolution, resolution, resolution, isolevel, V_out, F_out);

    if (V_out.rows() == 0 || F_out.rows() == 0) {
        std::cerr << "Marching Cubes failed: Empty mesh" << std::endl;
        return;
    }

    // calculate normal
    Eigen::MatrixXd N;
    igl::per_face_normals(V_out, F_out, N);
    std::cout << "Generated mesh: " << V_out.rows() << " vertices, " << F_out.rows() << " faces" << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Model generation completed in " << duration.count() << " ms" << std::endl;
    return;


}

double ModelGenerator::combinedSDF(Eigen::Vector3d & p)
{
    double total_void = 0.0;
    int gaus_num = Kernels.size();
    for (size_t i = 0; i < gaus_num; i++) {
        //std::cout << Kernels[i].gaussian_fun(p) << std::endl;
        total_void += Kernels[i].gaussian_fun(p);
    }
    
    return  gauss_combined - total_void;  // 当前使用：负的空洞总和 
}

void ModelGenerator::show_model()
{
    std::cout << "show libigl viewer" << std::endl;
    igl::opengl::glfw::Viewer viewer;

    viewer.data().set_mesh(V_ini, F_ini);
    viewer.data().show_lines = true;   // 不显示网格线
    //viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色

    int id2 = viewer.append_mesh();
    Eigen::MatrixXd V_shifted = V_out;
    V_shifted.rowwise() += Eigen::RowVector3d(150, 0.0, 0.0);  // 向右移动 1.5 个单位

    viewer.data(id2).set_mesh(V_shifted, F_out);
    viewer.data(id2).set_colors(Eigen::RowVector3d(0.8, 0.1, 0.1));

    // 添加辅助点 (高斯核的中心)，设置为红色
   // viewer.data().add_points(kernel_points, Eigen::RowVector3d(1, 0, 0));
    viewer.data().point_size = 10; // 让点更显眼

    viewer.launch();

}