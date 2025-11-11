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
		SDF_gaussian(idx) = combinedSDF(p, Kernels, gauss_combined);
        SDF_out(idx) = smooth_IntersecSDF(SDF_ini(idx), -SDF_gaussian(idx), smooth_t);
        //SDF_out(idx) = differenceSDF(SDF_ini(idx), SDF_gaussian(idx));
        //if (SDF_gaussian(idx) < isolevel)
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
    
    if (Kernels.size() != 2)  cout << "more kernel" << endl;
	//std::vector<GaussianKernel> fat_curve = generate_tube(Kernels[0], Kernels[1], 0.5, 0.5);
    Eigen::VectorXd SDF_gaussian_tubes(grid_num);
    for (int idx = 0; idx < grid_num; ++idx) {
        Eigen::Vector3d p = GV.row(idx);
        //SDF_gaussian_tubes(idx) = max(SDF_gaussian(idx), generate_tube(p, Kernels[0], Kernels[1], 0.5, 0.5));
        SDF_gaussian_tubes(idx) = smooth_UnionSDF(SDF_gaussian(idx), -generate_tube(p, Kernels[0], Kernels[1], gauss_combined, 0.5), smooth_t);
        //SDF_gaussian_tubes(idx) = generate_tube(p, Kernels[0], Kernels[1], gauss_combined, 0.5);
    }

    Eigen::MatrixXd V_t; //输出网格顶点
    Eigen::MatrixXi F_t; // 输出网格面片
    MarchingCubes(SDF_gaussian_tubes, GV, resolution, resolution, resolution, -0.1, V_t, F_t);
    view_model(V_t, F_t);

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

double ModelGenerator::combinedSDF(Eigen::Vector3d & p, std::vector<GaussianKernel> G_kernels, double C)
{
    double total_void = 0.0;
    int gaus_num = G_kernels.size();
    for (size_t i = 0; i < gaus_num; i++) {
        //std::cout <<"gaus_num: "<< gaus_num<<"  "<< i << std::endl;
        total_void += G_kernels[i].gaussian_fun(p);
    }
    //std::cout << "total_void: " << total_void << std::endl;
    return  C - total_void;  // 当前使用：负的空洞总和 
}

void ModelGenerator::show_model()
{
    view_two_models(V_ini, F_ini, V_out, F_out, Eigen::RowVector3d(1, 0.0, 0.0));

}


std::vector<Edge> ModelGenerator::pores_connection_mst(const std::vector<GaussianKernel>& gau)
{
    std::vector<Edge> mst_edges;

    int n = gau.size();
    if (n <= 1)
    {
		"Only one or no Gaussian kernels provided.";
        return mst_edges;
    }

    std::vector<bool> visited(n, false);
    std::priority_queue<Edge, std::vector<Edge>, CompareEdge> pq;

    // 从第一个节点开始
    visited[0] = true;

    // 初始化：加入与第一个节点相连的所有边
    for (int j = 1; j < n; ++j) {
        double dist = (gau[0].center - gau[j].center).norm();
        pq.push({ 0, j, dist });
    }

    // 逐步扩展生成树
    while (!pq.empty() && mst_edges.size() < n - 1) {
        Edge e = pq.top();
        pq.pop();

        if (visited[e.to]) continue; // 避免重复访问

        // 接受该边
        mst_edges.push_back(e);
        visited[e.to] = true;

        // 将新加入节点的边放入队列
        for (int j = 0; j < n; ++j) {
            if (!visited[j]) {
                double dist = (gau[e.to].center - gau[j].center).norm();
                pq.push({ e.to, j, dist });
            }
        }
    }
    return mst_edges;
}

GaussianKernel ModelGenerator::make_sphere_gaussian(const Eigen::Vector3d& center, double radius, double targetC)
{
    GaussianKernel g;
    g.center = center;
    // 经验：把 sigma 设为 radius / 2.0（可调），以保证曲线足够平滑且半径处接近 targetC
    double sigma = std::max(radius / 2.0, 1e-6);
    g.sigma = sigma;
    double exponent = -(radius * radius) / (2.0 * sigma * sigma);
    // prevent underflow
    exponent = std::max(exponent, -745.0); // exp(-745) ~ 5e-324
    double A = targetC / std::exp(exponent);
    // A 是我们这里定义的 amplitude，使得 amplitude * exp(-r^2/(2 sigma^2)) == targetC
    g.amplitude = A;
    return g;
}

/*
 generate_tube:
  - k1, k2: 两个 GaussianKernel，代表两个孔的核（使用 .center 与 .sigma 信息）
  - optional parameters could be added; for clarity we use some defaults:
     * nControl = 5 (Bernstein degree 4)
     * min_mid_radius_factor: 控制中间最小半径作为端点半径的比例（防止过细）
  - 结果：把生成的沿线若干 GaussianKernel push 到 fat_curve_kernels（全局）
*/
//std::vector<GaussianKernel>  ModelGenerator::generate_tube(const GaussianKernel& k1, const GaussianKernel& k2, double C, double mid_radius_factor ) // 中间最小半径相对端点半径的初值
//{
//    // degree n = 4 -> 5 control samples
//    const int n = 4;
//    const int sampleCount = n + 1;
//
//    Eigen::Vector3d p0 = k1.center;
//    Eigen::Vector3d p1 = k2.center;
//
//    // 端点“孔的内切球半径”的近似：根据 sigma 推断（可根据具体实现替换为更精确值）
//    // 我们认为：若 Gaussian 在距离 r0 处的值等于 C，则 r0 与 sigma 存在关系。
//    // 但此处如果 k1.sigma 存在，我们可以近似取端点半径 = 2 * k.sigma （经验值）
//    double r0 = std::max(2.0 * k1.sigma, 1e-4);
//    double r4 = std::max(2.0 * k2.sigma, 1e-4);
//
//    // 设置中间最小半径（按论文：由可制造性/孔大小约束设置）
//    double r_mid = std::max(mid_radius_factor * std::min(r0, r4), 1e-4);
//
//    // 五个控制半径 r_i (i=0..4): 两端为 r0,r4，中间由用户或经验给定（我们将中点 i=2 设为 r_mid，其余线性插值）
//    std::vector<double> radii(sampleCount, 0.0);
//    radii[0] = r0;
//    radii[4] = r4;
//    // set middle control radii: simple strategy: [r0, interp1, r_mid, interp2, r4]
//    radii[2] = r_mid;
//    radii[1] = 0.5 * (r0 + r_mid);
//    radii[3] = 0.5 * (r_mid + r4);
//
//    std::vector<GaussianKernel> fat_curve_kernels;
//    // For each sample t_i = i/(sampleCount-1), compute sample position and radius via Bernstein
//    for (int si = 0; si < sampleCount; ++si)
//    {
//        double t = double(si) / double(sampleCount - 1); // 0..1
//        // position = p0 * (1 - t) + p1 * t  (we use straight line; paper samples uniform along segment)
//        Eigen::Vector3d pos = p0 * (1.0 - t) + p1 * t;
//
//        // compute radius via Bernstein blending of control radii
//        double R = 0.0;
//        for (int i = 0; i <= n; ++i) {
//            double B = bernstein_basis(i, n, t);
//            R += B * radii[i];
//        }
//
//        // create a Gaussian kernel approximating the sphere of radius R located at pos,
//        // and choose amplitude so that at distance R the Gaussian value equals C (paper: fat-curve surface value = C).
//        GaussianKernel g = make_sphere_gaussian(pos, R, C);
//
//        // append to global container
//        fat_curve_kernels.push_back(g);
//    }
//    return fat_curve_kernels;
//    // note: 论文重建时把“fat-curve 的隐式函数”与原 Gaussian 合并（取最大值）并在体素上执行 marching cubes。
//    // 这里我们只负责生成沿线的高斯核近似（显式球族），后续合并与提取在你的体素/场表示中完成。
//}


double ModelGenerator::generate_tube(const Eigen::Vector3d& p, const GaussianKernel& k1, const GaussianKernel& k2, double iso_level_C, double mid_radius_factor) // 中间最小半径相对端点半径的初值
{
    // degree n = 4 -> 5 control samples
    const int n = 4;
    const int sampleCount = n + 1;

    Eigen::Vector3d c1 = k1.center;
    Eigen::Vector3d c2 = k2.center;


    Eigen::Vector3d line_dir = c2 - c1;
    double line_length = line_dir.squaredNorm();

    // 2. 计算点在连接线上的投影参数t [0,1]
    double t_proj = (p - c1).dot(line_dir) / line_length;

    // 如果投影点不在线段上，返回0（不在管道内）
    if (t_proj < 0.0 || t_proj > 1.0) {
        return 0.0;
    }
    // 将t限制在[0, 1]区间，这样空间中所有的点都会被映射到线段上最近的点
    t_proj = std::max(0.0, std::min(1.0, t_proj));
    Eigen::Vector3d p_proj = c1 + t_proj * line_dir; // p在轴线上的投影点

    double r0 = k1.sigma * std::sqrt(2.0 * std::log(k1.amplitude / iso_level_C));
    double r4 = k2.sigma * std::sqrt(2.0 * std::log(k2.amplitude / iso_level_C));

    // 设置中间最小半径（按论文：由可制造性/孔大小约束设置）
    double r2 = std::max(mid_radius_factor * std::min(r0, r4), 1e-4);

    double r1 = (r0 + r2) / 2.0;
    double r3 = (r2 + r4) / 2.0;

    std::vector<double> control_radii = { r0, r1, r2, r3, r4 };

    // 使用Bernstein多项式计算R(t)
    double R_t = 0.0;
    for (int i = 0; i <= 4; ++i) {
        R_t += bernstein_basis(i, 4, t_proj) * control_radii[i];
    }

    // 如果目标半径非常小，则认为没有通道贡献
    //if (R_t < 1e-6) {
    //    return 0.0;
    //}
    // 3. 计算点p的隐式函数值
// 核心思想是：函数值应该在通道表面上为C，在轴线上最大，向外围衰减。
    double dist_to_axis = (p - p_proj).norm();
    double energy = (dist_to_axis * dist_to_axis) / (R_t * R_t);
    //return 1.0 - energy;
    return iso_level_C * std::exp(1.0 - energy);



    //if (dist_to_axis <= R_t) {
    //    // 在管道内部，使用平滑的隐式函数
    //    // 这里使用二次函数，在边界处值为0，在轴线处值为1
    //    double normalized_dist = dist_to_axis / R_t;
    //    return 1.0 - normalized_dist * normalized_dist;
    //}
    //else {
    //    return 0.0;  // 在管道外部
    //}
}