#include "modelGen.h"

//GaussianKernel::GaussianKernel(Eigen::Vector3d center_, double sigma_, double amplitude_, double center_value_)
//{
//    center = center_;
//    sigma = sigma_;         // 核的大小/影响力范围 (高斯函数的标准差)
//    amplitude = amplitude_;
//    center_value = center_value_;
//    on_surface = is_on_surface();  //true: on surface
//}
//
//double GaussianKernel::gaussian_fun(const Eigen::Vector3d& p)
//{
//    // ---------- 原始高斯值 ----------
//    double d2 = (p - center).squaredNorm();
//    double G = amplitude * std::exp(-d2 / (2.0 * sigma * sigma));
//    return G;
//}
//
//bool GaussianKernel::is_on_surface() const
//{
//    double ratio = Gauss_level / amplitude;
//    double d_value = sqrt(-2.0 * sigma * sigma * std::log(ratio));
//    //cout << d_value + 1e-3 << "   " << abs(center_value) << endl;
//    return d_value + 1e-3 > abs(center_value);
//}

GaussianKernel::GaussianKernel(
    const Eigen::Vector3d& center_,
    double sigma_parallel_,
    double sigma_perp_,
    const Eigen::Matrix3d& rotation_,
    double amplitude_,
    double center_value_
)
{
    center = center_;
    sigma_parallel = sigma_parallel_;
    sigma_perp = sigma_perp_;
    R = rotation_;
    amplitude = amplitude_;
    center_value = center_value_;
    on_surface = is_on_surface();

    // 局部坐标下的 Σ^{-1}
    Eigen::Matrix3d Dinv = Eigen::Matrix3d::Zero();
    Dinv(0, 0) = 1.0 / (sigma_perp * sigma_perp);
    Dinv(1, 1) = 1.0 / (sigma_perp * sigma_perp);
    Dinv(2, 2) = 1.0 / (sigma_parallel * sigma_parallel);
    //R = Eigen::Matrix3d::Identity(); //暂时只沿建造方向拉长

    // 世界坐标下的 Σ^{-1} = R * D^{-1} * R^T
    invSigma = R * Dinv * R.transpose();
}

double GaussianKernel::gaussian_fun(const Eigen::Vector3d& p) const
{
    Eigen::Vector3d d = p - center;
    double quad = d.transpose() * invSigma * d;
    return amplitude * std::exp(-0.5 * quad);
}

bool GaussianKernel::is_on_surface() const
{
    double ratio = Gauss_level / amplitude;
    double d_factor = std::sqrt(-2.0 * std::log(ratio));
    // 对各向异性椭球，取最大轴作为判断
    double d_max = d_factor * std::max(sigma_parallel, sigma_perp);

    return d_max + 1e-3 > abs(center_value);
}

ModelGenerator::ModelGenerator(std::string input_file, int pores)
{
    if (!igl::read_triangle_mesh(input_file, V_ini, F_ini)) {
        std::cerr << "Error: Could not load model A." << std::endl;
        return;
    }
    pore_num = pores;
}

void ModelGenerator::model_porous_structure()
{
    //std::cout << "Model A loaded successfully." << std::endl;
    scale_factor = Mesh2SDF(V_ini, F_ini, GV, SDF_ini, bb_min, bb_max);
    //Vector3d point = GV.row(10010);
    //cout << "point index: " << find_nearest_grid(point) << endl;
    generateGaussianSDF();


    /*VoxelGrid grids = SDFtoVoxel(SDF_out, bb_min, bb_max, resolution, resolution, resolution);
    SupportCheckResult scr = checkSupportVoxel(grids, 0.5);*/


    //// 保存为NPY格式
    //std::string npy_filename = "D:/VSprojects/TaihuStone/model/npy/voxelized_model_out.npy";
    //saveVoxelGridAsNPY(voxel_grid, Resolution, npy_filename);
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
    std::vector<double> pore_sdfs;
    std::vector<int> inside_indices;
    sample_interior_points(pore_centers, pore_sdfs, inside_indices, pores, gen);
    double model_count = inside_indices.size();
    //for (int i = 0; i < pore_centers.size(); i++)   std::cout << "i: " << i << "  " << pore_centers[i] << std::endl;

    generate_gaussians(pore_centers, pore_sdfs, gen);
    
	int kernel_num = Kernels.size();
    int degree_limit = (kernel_num - 1) / max(1, (int)(kernel_num - surface_kernels.size()));

    degree_limit = max(degree_limit, Min_degree);
    Tube_edges = pores_connection_mst(Kernels, degree_limit);

	std::vector<pair<int, int>> edge_connections;
	for (auto e : Tube_edges) edge_connections.push_back(make_pair(e.from, e.to));
    vis_Kernels_Tubes(pore_centers, edge_connections);
	// 构建邻接表
    Adj_list = construct_adj_list(Tube_edges, kernel_num);
    Unused_adj_list = get_unused_edge_adj(Adj_list, Adj_dis_thres);
    vector<int> leafs_index = all_leafs_mst(Tube_edges);
    vector<int> inner_leafs = check_inner_leafs(leafs_index);
    
	//---------calculate total translucency score------------------------------
    std::cout << "--------------------3. Calculating total translucency of mst --------------------" << endl;
    finalTranslucency = cal_total_translucency(Kernels, Adj_list);
	int ori_edge_num = Tube_edges.size();
    std::cout << "Before optimization, total score: " << finalTranslucency << " with "<< ori_edge_num <<" edges."<<endl;

	//--------------optimize connection edges------------------------------

    //vector<int> edge_improtance = cal_edge_usage(Paths, true);
	// 输出每条边的重要性分数
	//for (auto ei : edge_improtance) 
 //       cout << "edge_importance: " << ei << endl;
    std::cout << "--------------------4. Optimizing the connection trees --------------------" << endl;
    int opt_times_once = 5;
	int edge_max = Tube_edges.size()* 1.5;
    int iter_times = 10;
    double delta_score_t = 1000.0;
    double new_trans_score = finalTranslucency;
    
    if (optimize_debug)
    {
        //optimize_mst2(opt_times_once, edge_max, false,  NO_DEBUG);
		int iter_count = 0;
        while (iter_count< iter_times && delta_score_t > 1e-2)
        {
            std::cout << "-------------  The " << iter_count++<< " iterations start, the edge limition is:"<< edge_max<<"... ------------ - "<< endl;
            optimize_mst(opt_times_once, edge_max, NO_DEBUG);
            double last_trans_score = new_trans_score;
            new_trans_score = cal_total_translucency(Kernels, Adj_list);
            std::cout << "======================================After optimization " << iter_count << ", total score increases from " << last_trans_score << " to " << new_trans_score << " with edges to " << Tube_edges.size() << "======================================" << endl;
            delta_score_t = new_trans_score - last_trans_score;
            edge_max += edge_max * 0.1;
        }
    }
    finalTranslucency = new_trans_score; 
    std::vector<pair<int, int>> edge_conn_new;
    for (auto e : Tube_edges) edge_conn_new.push_back(make_pair(e.from, e.to));
    vis_Kernels_Tubes(pore_centers, edge_conn_new);
    

    std::cout << "--------------------5. Generate tubes between kernels based on optimized mst --------------------" << endl;
    //-----------------generate tubes------------------------------------------
    double void_count = generate_mst_tubes(grid_num, resolution, Isolevel, Gauss_level, SmoothT);
    //std::string vti_filename = "D:/VSprojects/TaihuStone/result/vti/" + input_file + "_voxelized_model_" + std::to_string(resolution) + "x" + std::to_string(resolution) + "x" + std::to_string(resolution) + ".vti";
    //saveSDFtoVTI(vti_filename, SDF_out, resolution, resolution, resolution);

    std::string raw_filename = "D:/VSprojects/TaihuStone/result/raw/" + input_file + "_voxelized_model_" + std::to_string(resolution) + "x" + std::to_string(resolution) + "x" + std::to_string(resolution) + ".vol";
    saveSDFtoVOI(raw_filename, SDF_out, resolution, resolution, resolution);

    VoxelGrid grids = SDFtoVoxel(SDF_out, bb_min, bb_max, resolution, resolution, resolution);
    //SupportCheckResult scr = checkSupportVoxel(grids, 0.5);

    std::string npy_filename = "D:/VSprojects/TaihuStone/model/npy/"+ input_file+ "_voxelized_model_" + std::to_string(resolution) + "x" + std::to_string(resolution) + "x" + std::to_string(resolution) + ".npy";
    saveVoxelGridAsNPY(grids.rho, resolution, npy_filename);

    // 计算孔隙率
    double porosity = 1.0 - void_count / model_count;
    std::cout << "Porosity: " << porosity * 100 << "%" << "    --------:" << void_count << "   " << model_count << std::endl;

    // 存储最终的孔隙率
    finalPorosity = porosity;

    //out file
    std::string filename = "result/porous_" + input_file;
    if (optimize_debug)
        filename += "_opt";
    filename = filename + "_" + to_string(PoresNum) + "_" + to_string(Resolution) + "_" + to_string_pre(surface_ratio) + "_" +
        to_string_pre(Trans_thres) + "_" + to_string_pre(Weights[0]) + "_" + to_string_pre(KT_weights[0]) + "_" +
        to_string_pre(sigma_min, 3) + "_" + to_string_pre(sigma_max, 3) + "_" + to_string_pre(finalTranslucency, 3) + "_"+ to_string_pre(finalPorosity*100.0)+"%.stl";
    saveMesh(filename, V_out, F_out);
	// ------------check single component--------------
    //int comp = single_component(V_g, F_g);
  

    return;

}

void ModelGenerator::sample_interior_points(std::vector<Eigen::Vector3d>& pore_centers, std::vector<double>& pore_sdfs, std::vector<int>& inside_indices, int pores, std::mt19937& gen)
{
    Eigen::VectorXd SDF = this->SDF_ini;
    Eigen::MatrixXd GV = this->GV;
    int grid_num = SDF.size();
	double margin = 0.03;
    double out_margin = 0.3 * margin;
    std::vector<int> surface_indices;

    // search inside points
    for (int idx = 0; idx < grid_num; ++idx) {
        if (SDF(idx) < Isolevel) {
            inside_indices.push_back(idx);
        }
        if (SDF(idx) < out_margin && SDF(idx) > -margin) {   //找到边界区域
            surface_indices.push_back(idx);
        }
    }
    double model_count = inside_indices.size();
	double surface_count = surface_indices.size();
	cout << "Total grid num: " << grid_num << ", inside grid num: " << model_count << ", surface grid num: " << surface_count <<endl;
    if (inside_indices.empty()) {
        std::cerr << "Warning: no legal points inside!" << std::endl;
        return;
    }

    // 在整个形状内部采样
    std::uniform_int_distribution<int> index_dist(0, inside_indices.size() - 1);
    std::uniform_int_distribution<int> surface_dist(0, surface_indices.size() - 1);

    int attempts = 0;
    int max_attempts = pores * 50;
    double suf_ratio = surface_ratio;
	int surface_sam = pores * suf_ratio;
    int surface_p = 0;
    int all_sam_num = pore_centers.size();
    Eigen::Vector3d min_pt = GV.colwise().minCoeff();
    Eigen::Vector3d max_pt = GV.colwise().maxCoeff();
    Eigen::Vector3d box_size = max_pt - min_pt;
    double volume = box_size.x() * box_size.y() * box_size.z();
    double safe_unit = std::cbrt(volume / pores);
    safe_distance = Safe_distance_ratio * safe_unit;
    int base_layer = find_nearest_grid(bb_min) / 10000;

    while (all_sam_num < pores && attempts < max_attempts) {
        attempts++;
        if (attempts > max_attempts - 2) //逼近界限，安全距离缩减，提升20次上限
        {
            Safe_distance_ratio = Safe_distance_ratio - 0.01;
            safe_distance = Safe_distance_ratio * safe_unit; 
            max_attempts = max_attempts + 50;
        }
        // 随机选择一个内部点，保存其sdf值
        int chosen_idx = -1;
        if(surface_p< surface_sam)
            chosen_idx = surface_indices[surface_dist(gen)];
        else{
            chosen_idx = inside_indices[index_dist(gen)];
            if (SDF(chosen_idx) < out_margin && SDF(chosen_idx) > -margin)
                    continue;
        }
        //cout << chosen_idx << "    " << chosen_idx / (Resolution * Resolution) << "   " << base_layer + 15 << endl;
        //if (chosen_idx / (Resolution * Resolution) < base_layer + 10)
        //{
        //    cout << chosen_idx << "    " << chosen_idx / (Resolution * Resolution) << "   " << base_layer + 10 << endl;
        //    continue;
        //}
        Eigen::Vector3d candidate_center = GV.row(chosen_idx).transpose();
        // 检查与已有空洞中心的最小距离
        bool valid = true;
        for (const auto& existing_center : pore_centers) {
            if ((candidate_center - existing_center).squaredNorm() < safe_distance * safe_distance) {
                valid = false;
                break;
            }
        }

        if (valid) {
            pore_centers.push_back(candidate_center);
            pore_sdfs.push_back(SDF(chosen_idx));
            all_sam_num++;
            if (SDF(chosen_idx) < out_margin && SDF(chosen_idx) > -margin)
                surface_p++;
            if(debug_show)
                cout << "pore_sdfs: " << SDF(chosen_idx) << endl;
        }
    }

    std::cout << "Generate " << pore_centers.size() << " kernels   " << all_sam_num<<" , include "<< surface_p <<" sp with safe_dis: "<< Safe_distance_ratio << std::endl;
}

void ModelGenerator::generate_gaussians(std::vector<Eigen::Vector3d> pore_centers, std::vector<double> pore_sdfs, std::mt19937& gen)
{
    //Kernels.clear();
    //surface_kernels.clear();
    //bool dynamic_change_para = true;
    //if(dynamic_change_para)
    //{
    //    double bbx_v = abs((bb_max.x() - bb_min.x()) * (bb_max.y() - bb_min.y()) * (bb_max.z() - bb_min.z()));
    //    double w4max = 7.0, w4min = 42.0;
    //    sigma_min = sigma_value(bbx_v, pore_num, w4min, Gauss_level);
    //    sigma_max = sigma_value(bbx_v, pore_num, w4max, Gauss_level);
    //    cout << "sigma_min: " << sigma_min << "   sigma_max: " << sigma_max << endl;
    //}
    //
    //std::uniform_real_distribution<double> dist_amp(amplitude_min, amplitude_max);
    //std::uniform_real_distribution<double> dist_sigma(sigma_min, sigma_max);

    //// 为每个空洞中心生成随机参数
    //int pore_size = pore_centers.size();
    //for (size_t i = 0; i < pore_size; ++i) {
    //    double sigma_val = dist_sigma(gen);
    //    double amplitude_val = dist_amp(gen);
    //    //std::cout << "i: " << i << "  " << sigma_val<<"  "<<amplitude_val << std::endl;
    //    GaussianKernel kernel(pore_centers[i], sigma_val, amplitude_val, pore_sdfs[i]);
    //    Kernels.push_back(kernel);
    //    if (kernel.on_surface) surface_kernels.push_back(i);
    //    /*pore_amplitudes.push_back(dist_amp(gen));
    //    pore_sigmas.push_back(dist_sigma(gen));*/
    //}

    //std::cout << "--------------------1. Sample and generate gaussian kernels--------------------" << endl<<
    //    "Combine " << Kernels.size() << " Gaussian fileds with " << surface_kernels.size() << " kernels on the surface" << std::endl;
    //if(standard_show)
    //{
    //    std::cout << "surface index:  ";
    //    for (auto i : surface_kernels) std::cout << i << "  ";
    //    std::cout << endl;
    //}
    Kernels.clear();
    surface_kernels.clear();
    bool dynamic_change_para = true;
    if (dynamic_change_para)
    {
        double bbx_v = abs((bb_max.x() - bb_min.x()) * (bb_max.y() - bb_min.y()) * (bb_max.z() - bb_min.z()));
        double w4max = 7.0, w4min = 42.0;
        sigma_min = sigma_value(bbx_v, pore_num, w4min, Gauss_level);
        sigma_max = sigma_value(bbx_v, pore_num, w4max, Gauss_level);
        cout << "sigma_min: " << sigma_min << "   sigma_max: " << sigma_max << endl;
    }

    std::uniform_real_distribution<double> dist_amp(amplitude_min, amplitude_max);
    std::uniform_real_distribution<double> dist_sigma(sigma_min, sigma_max);

    // 为每个空洞中心生成随机参数
    int pore_size = pore_centers.size();
    for (size_t i = 0; i < pore_size; ++i) {
        double sigma_val = dist_sigma(gen);
        double amplitude_val = dist_amp(gen);
        //std::cout << "i: " << i << "  " << sigma_val<<"  "<<amplitude_val << std::endl;
        //GaussianKernel kernel(pore_centers[i], sigma_val * 1.5, sigma_val, Eigen::Matrix3d::Identity(), amplitude_val, pore_sdfs[i]);
        GaussianKernel kernel(pore_centers[i], sigma_val, sigma_val, Eigen::Matrix3d::Identity(), amplitude_val, pore_sdfs[i]);
        Kernels.push_back(kernel);
        if (kernel.on_surface) surface_kernels.push_back(i);
        /*pore_amplitudes.push_back(dist_amp(gen));
        pore_sigmas.push_back(dist_sigma(gen));*/
    }

    std::cout << "--------------------1. Sample and generate gaussian kernels--------------------" << endl <<
        "Combine " << Kernels.size() << " Gaussian fileds with " << surface_kernels.size() << " kernels on the surface" << std::endl;
    if (standard_show)
    {
        std::cout << "surface index:  ";
        for (auto i : surface_kernels) std::cout << i << "  ";
        std::cout << endl;
    }
}


double ModelGenerator::combinedSDF(Eigen::Vector3d & p, std::vector<GaussianKernel> G_kernels, double C)
{
    double total_void = 0.0;
    int gaus_num = G_kernels.size();
    for (size_t i = 0; i < gaus_num; i++) {
        //std::cout <<"gaus_num: "<< gaus_num<<"  "<< i << std::endl;
        total_void += G_kernels[i].gaussian_fun(p);
    }

    //double noise_val = myPerlin.noise(p.x() * 0.5, p.y() * 0.5, p.z() * 0.5);
    //double noise_weight = 0.05; // 权重

    //// 注意：SDF通常定义为 (IsoLevel - FieldValue)，所以加在 FieldValue 上
    //total_void = total_void *(1+ noise_weight * noise_val) ;
    //std::cout << "total_void: " << total_void << std::endl;
    return  C - total_void;  // 当前使用：负的空洞总和 
}

void ModelGenerator::show_model()
{
    std::cout << "show libigl viewer" << std::endl;
    igl::opengl::glfw::Viewer viewer;
    viewer.core().background_color.setConstant(1.0f); // White background
    Eigen::RowVector3d shift(1, 0, 0);

    viewer.data().set_mesh(V_out, F_out);
    //viewer.data().show_lines = false;   // 不显示网格线
    viewer.data().show_lines = true;   // 不显示网格线
    viewer.data().show_faces = false;   // 不显示三角面
    viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.8, 0.8));
    viewer.data().shininess = 500.0;
    //viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色

    int id2 = viewer.append_mesh();
    Eigen::MatrixXd V_shifted = V_out;
    V_shifted.rowwise() -= shift;
    viewer.data(id2).show_faces = false;   // 不显示三角面
    viewer.data(id2).set_mesh(V_shifted, F_out);
    viewer.data(id2).shininess = 500.0;
    viewer.data(id2).set_colors(Eigen::RowVector3d(0.8, 0.8, 0.8));

    int id3 = viewer.append_mesh();
    viewer.data(id3).set_mesh(V_t, F_t);
    //viewer.data(id3).set_colors(Eigen::RowVector3d(0.4, 0.4, 0.2));
    viewer.data(id3).show_lines = false;   // 不显示网格线
    viewer.data(id3).shininess = 500.0;
    Eigen::MatrixXd V_shifted3 = V_t;
    V_shifted3.rowwise() += shift;  // 向右移动 1 个单位

    int id4 = viewer.append_mesh();
    viewer.data(id4).set_mesh(V_shifted3, F_t);
    viewer.data(id4).show_lines = false;   // 不显示网格线
    viewer.data(id4).shininess = 500.0;
    // 添加辅助点 (高斯核的中心)，设置为红色

    viewer.data().point_size = 10; // 让点更显眼
    viewer.launch();
    //
    //view_two_models(V_ini, F_ini, V_out, F_out, Eigen::RowVector3d(1, 0.0, 0.0));

}


std::vector<Edge> ModelGenerator::pores_connection_mst(const std::vector<GaussianKernel>& gau, int Dmax)
{
    cout << "--------------------2. Constructing connection MST with most degree Dmax: " << Dmax << "--------------------"<<endl;
    std::vector<Edge> mst_edges;

    int n = gau.size();
    if (n <= 1)
    {
        std::cout << "Only one or no Gaussian kernels provided." << std::endl;
        return mst_edges;
    }

    std::vector<bool> visited(n, false);
    std::vector<int> degree(n, 0);   // <-- 新增：记录节点度数
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;

    // 从第一个节点开始
    visited[0] = true;

    // 初始化：加入与第一个节点相连的所有边
    for (int j = 1; j < n; ++j) {
        double dist = distance(gau[0].center, gau[j].center);
        double dist_w = dist * calculate_edge_weight(gau[0], gau[j]);
        pq.push({ 0, j, dist, dist_w });
    }

    // 逐步扩展生成树
    while (!pq.empty() && mst_edges.size() < n - 1) {
        Edge e = pq.top();
        pq.pop();

        if (visited[e.to]) continue; // 避免重复访问

        // 度数限制检查
        if (degree[e.from] >= Dmax || degree[e.to] >= Dmax)
            continue;
        // 接受该边
        mst_edges.push_back(e);
        visited[e.to] = true;

        degree[e.from]++;
        degree[e.to]++;

        // 将新加入节点的边放入队列
        for (int j = 0; j < n; ++j) {
            if (!visited[j]) {
                double dist = distance(gau[e.to].center, gau[j].center);
                double dist_w = dist * calculate_edge_weight(gau[e.to], gau[j]);
                pq.push({ e.to, j, dist, dist_w });
            }
        }
    }
    if(standard_show)
        for (Edge e : mst_edges)
            std::cout << "Edge from " << e.from << " to " << e.to << " with weight " << e.weight << " with length " << e.length << std::endl;
    return mst_edges;
}

std::vector<std::vector<int>> ModelGenerator::construct_adj_list(std::vector<Edge> edges_list, int kernel_num)
{
    std::vector<std::vector<int>> adj_list(kernel_num);
    for (const auto& edge : edges_list) {
        // 由于最小生成树是无向图，一条边代表双向连接, 需要将 `to` 添加到 `from` 的邻居列表，同时将 `from` 添加到 `to` 的邻居列表。
        if (edge.from < kernel_num && edge.to < kernel_num) {
            adj_list[edge.from].push_back(edge.to);
            adj_list[edge.to].push_back(edge.from);
        }
    }
	return adj_list;
}

std::vector<std::vector<int>> ModelGenerator::get_unused_edge_adj(AdjacencyList Adj_list, double dis_thres)
{
    int n = Kernels.size();

    // 初始化结果邻接表
    std::vector<std::vector<int>> unused_adj(n);

    // 预计算阈值的平方，避免循环中开根号
    double threshold_sq = dis_thres * dis_thres;

    // 遍历所有唯一的点对 (i, j)
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            // 1. 距离筛选 (最快，先做)
            double dist_sq = squared_distance(Kernels[i].center, Kernels[j].center);

            if (dist_sq > threshold_sq) {
                continue; // 距离过大，剔除
            }
            // MST 存在性检查
            bool is_in_mst = false;
            for (int neighbor : Adj_list[i])
            {
                if (neighbor == j) {
                    is_in_mst = true;
                    break;
                }
            }
            if (is_in_mst) {
                continue; // 边已被 MST 使用，剔除
            }
            // 3. 添加到结果 (双向)
            unused_adj[i].push_back(j);
            unused_adj[j].push_back(i);
        }
    }
    return unused_adj;
}

double ModelGenerator::cal_path_graph_length(std::vector<int> path_)
{
    double length = 0;
    for (size_t i = 0; i < path_.size() - 1; ++i) {
        int u = path_[i];
        int v = path_[i + 1];
        length += length_path(u, v);
	}
    return length;
}


std::vector<int> ModelGenerator::find_path_in_tree(int p1, int p2,  int num_nodes, AdjacencyList adj)   //BFS
{
    // --- 步骤 1: 将边列表转换为邻接表 ，其中 adj_list[i] 将存储所有与节点 i 直接相连的节点的ID。
    if (num_nodes <= 0) {
        return {}; // 如果没有节点，则不可能有路径
    }
    if (p1 == p2) {
        // 处理起点和终点是同一个节点的特殊情况
        return {p1}; // 如果节点重复，则输出p1
    }
    std::vector<std::vector<int>> adj_list = adj;

    // --- 步骤 2: 使用广度优先搜索 (BFS) 查找路径 , BFS 是在树或图中查找最短路径（按边数计）的经典算法。
    std::queue<int> q; // q: 一个队列，用于存放待访问的节点
    q.push(p1);

    // parent: 一个数组，用于在找到路径后进行回溯。parent[i] 存储的是在路径上节点 i 的前一个节点。
    std::vector<int> parent(num_nodes, -1);

    // visited: 一个布尔数组，用于标记节点是否已被访问，防止在图中走回头路或陷入循环。
    std::vector<bool> visited(num_nodes, false);
    visited[p1] = true;

    bool path_found = false;
    while (!q.empty()) {
        int u = q.front();
        q.pop();

        if (u == p2) {
            path_found = true;
            break; // 找到了终点，可以提前结束搜索
        }

        // 遍历当前节点 u 的所有邻居
        for (int v : adj_list[u]) {
            if (!visited[v]) {
                visited[v] = true; // 标记为已访问
                parent[v] = u;    // 记录父节点，用于后续路径重构
                q.push(v);        // 将邻居节点加入待访问队列
            }
        }
    }

    // --- 步骤 3: 路径重构:如果找到了路径（即BFS到达了终点）,从终点开始, 利用 `parent` 数组反向回溯，直到回到起点。
    std::vector<int> path;
    if (path_found) {
        int current = p2;
        while (current != -1) {
            path.push_back(current);
            current = parent[current]; // 移动到路径上的前一个节点
        }
        // 因为我们是-从终点回溯到起点，所以得到的路径是反的，需要反转一次。
        std::reverse(path.begin(), path.end());
    }

    return path;
}


std::vector<int> ModelGenerator::find_specified_path(int p_index, int s1, int s2, AdjacencyList adj)
{
    if (s1 == s2 || s1 < 0 || s2 < 0)
    {
        cout << "Warnning: cannot find an illegal path!" << endl;
        return {};
    }
    int kernel_num = Kernels.size();

    std::vector<int> path1 = find_path_in_tree(p_index, s1, kernel_num, adj);
    std::vector<int>  path2 = find_path_in_tree(p_index, s2, kernel_num, adj);

    // 合并路径，避免重复包含 p_index
    std::vector<int> full_path(path1.rbegin(), path1.rend());
    full_path.insert(full_path.end(), path2.begin() + 1, path2.end()); // 从 path2 的倒数第二个元素开始添加
    bool debug_show = false;
    if (debug_show)
    {
        cout << "Kernel " << p_index << "'s full path steps : " << full_path.size() << endl;
        for (auto p : full_path) cout << p << " ";
        cout << endl;
    }
    return full_path;
}

double ModelGenerator::length_graph_path(int p1, int p2, AdjacencyList adj)
{
    std::vector<Edge> mst = Tube_edges;
    if (p1 == p2) return 0.0;
    if (mst.empty()) return -1.0; // 错误：没有树

    std::vector<int> path_ = find_path_in_tree(p1, p2, Kernels.size(), adj);

	double length_graph = cal_path_graph_length(path_);
    return length_graph; // 如果图不连通（理论上MST应该是连通的），返回-1
}

double ModelGenerator::length_path(int p1, int p2)
{
    return distance(Kernels[p1].center, Kernels[p2].center);
}

int  ModelGenerator::find_edge_by_nodes(int from_node, int to_node, const std::vector<Edge> edge_list)
{
    double edge_index = -1;
    for (int i = 0; i < edge_list.size(); i++)
    {
        if ((edge_list[i].from == from_node && edge_list[i].to == to_node) ||
            (edge_list[i].from == to_node && edge_list[i].to == from_node))
        {
            edge_index = i;
        }
    }
    return edge_index;
}

std::vector<int> ModelGenerator::all_leafs_mst(std::vector<Edge>& mst_tree) 
{
    // 使用哈希表存储每个节点的度数
    // Key: 节点ID, Value: 度数
    std::unordered_map<int, int> node_degrees;

    // 1. 遍历所有边，统计度数
    for (const auto& edge : mst_tree) {
        node_degrees[edge.from]++;
        node_degrees[edge.to]++;
    }

    std::vector<int> leaves;

    // 2. 找出度数为 1 的节点
    for (const auto& pair : node_degrees) {
        if (pair.second == 1) {
            leaves.push_back(pair.first);
        }
    }

    // 3. 对结果进行排序（可选，但通常为了结果确定性建议排序）
    std::sort(leaves.begin(), leaves.end());

    return leaves;
}

bool ModelGenerator::find_edge_in_path(Edge cand_edge, vector<int> path)
{
    if (path.empty()) 
    {
        cout << "Empty path!" << endl;
        return true;  //如果路径为空，默认为在
    }
    //show_path(path);
    bool in = false;
    int sd = cand_edge.from;
	int ed = cand_edge.to;
    for (int ii = 0; ii < path.size()-1; ii++)
    {
        if (path[ii] == sd || path[ii] == ed)
        {
            if (path[ii + 1] == sd || path[ii + 1] == ed)
            {
                in = true;
				break;
            }
		}
    }
    return in;
}


double ModelGenerator::calculate_edge_weight(GaussianKernel k1, GaussianKernel k2)
{
    vector<double> weights = Weights;  // 分类系数权重：边界-内部，内部-内部，边界-边界
    // 1. 确定分类系数 C(i, j)
    double connect_weight;
    if (k1.on_surface != k2.on_surface) {
        // 一个是内部，一个是边界
        connect_weight = weights[0];
    }
    else if (!k1.on_surface) {
        // 两个都是内部
        connect_weight = weights[1];
    }
    else { // node1.type == NodeType::Boundary
        // 两个都是边界
        connect_weight = weights[2];
    }

    double overlap_weight = 1.0;
    // 2. 计算基础欧氏距离 d(i, j)


    // 3. 计算重叠度调节因子 O(i, j)
    //double O = calculate_overlap_factor(node1, node2, k_o);

    // 4. 计算最终权重 W(i, j) = d * C * O
    return connect_weight * overlap_weight;
}


double ModelGenerator::calculate_path_translucency(std::vector<int>& path, bool show_debug)
{
    std::vector<GaussianKernel> all_nodes = Kernels;
    int psize = path.size();
    double angle_product = 1.0;         // ∈ (0,1]
    int count_inner = 0;
    int count_surface = 0;
	int count_surface_line = 0;
	double alpha = 0.5;
    double L0 = 5.0;  //通透的标准长度
    double beta = 8.0;
    double mu = 0.0;
    double w_angle = KT_weights[0];
    double w_length = KT_weights[1];
    double w_location = KT_weights[2];
    double w_direction = KT_weights[3];
    std::vector<Eigen::Vector3d> path_points;
    
    if (psize < 2) {
        cout << "Warnning: illegal path!" << endl;
        return 0.0; // 单个点
    }
    //for (auto p : path) cout << p << "   ";
    //cout << endl;
    for (auto p : path) path_points.push_back(all_nodes[p].center);

	//get basic information
    for (size_t i = 1; i < psize - 1; ++i)
    {
        Vector3d prev = path_points[i - 1];
        Vector3d curr = path_points[i];
        Vector3d next = path_points[i + 1];
        if (!all_nodes[path[i]].on_surface) 
            count_inner++;
        double thres = min(all_nodes[path[i - 1]].center_value, all_nodes[path[i]].center_value);
        //line_cross_surface返回1.0，贯通； <1.0，属于表面，计数
        if (line_cross_surface(prev, curr, thres) < 0.999)
            count_surface_line++;
        if (i == psize - 2)
        {
            thres = min(all_nodes[path[i]].center_value, all_nodes[path[i + 1]].center_value);
            if (line_cross_surface(curr, next, thres) < 0.999)
                count_surface_line++;
        }
        double angle_deg = abs_angle(prev - curr, next - curr) / 180.0;
        angle_product *= angle_deg;
    }

    // ---------- 1. Angle term----------
    double T_angle;
    const double angle_floor = 0.6; // η_angle

    if (psize == 2) {
        T_angle = angle_floor;
    }
    else {
        T_angle = std::pow(angle_product, alpha);
    }

    // ---------- 2. Length term ----------
    double T_length = 1.0 - std::exp(-double(psize-1) / L0);


    // ---------- 3. Location term ----------
    double r_inner = (psize == 2)? 0: (double(count_inner) / double(psize - 2));
    double r_surface = (psize == 2) ? 0 : (double(count_surface_line) / double(psize - 1));

    double T_location = 1.0 / (1.0 + std::exp(-beta * (r_inner - r_surface - mu)));
	//cout << "T_location: "<<T_location << " = " <<"1/(1+e^-" <<beta << "  * (" << r_inner << " - " << r_surface << "  - " << mu <<"))" <<endl;
    // ---------- 4. Direction term ----------
    Eigen::Vector3d z(0, 0, 1);
    Eigen::Vector3d dir = computePrincipalDirection(path_points);
    double S_horiz = 1.0 - std::abs(dir.dot(z));
    //cout << "dir: " << dir <<"   S_horiz:  "<< S_horiz<< endl;
    double translucency_score = w_angle * T_angle + w_length * T_length + w_location * T_location + w_direction * S_horiz;
    //cout << translucency_score << "  " << T_angle << "   " << T_length << "  " << T_location << "   " << S_horiz << endl;
    return translucency_score;
}

double ModelGenerator::calculate_path_translucency2(std::vector<int>& path, bool show_debug)
{
    double translucency_score = 1.0;
    int psize = path.size();
    int sam_num = 3;
    if (psize < 2) {
        cout << "Warnning: illegal path!" << endl;
        return 0.0; // 单个点
    }
    if (psize == 2)
    {
        //是否横贯模型，如果是返回1，否则返回距离
        int start_ = path[0];
        int end_ = path[1];
        double thres = min(Kernels[start_].center_value, Kernels[end_].center_value);
        translucency_score = line_cross_surface(Kernels[start_].center, Kernels[end_].center, thres, sam_num);
    }
    else
    {
        double angle_product = 1.0;
        std::vector<GaussianKernel> all_nodes = Kernels;
        int count_inner = 0;
        int count_surface_line = 0;
        vector<double> angle_degrees;
        for (size_t i = 1; i < psize - 1; ++i)
        {
            Vector3d prev = all_nodes[path[i - 1]].center;
            Vector3d curr = all_nodes[path[i]].center;
            Vector3d next = all_nodes[path[i + 1]].center;
            if (!all_nodes[path[i]].on_surface) count_inner++;

            double thres = min(all_nodes[path[i - 1]].center_value, all_nodes[path[i]].center_value);
            //line_cross_surface返回1.0，贯通； <1.0，属于表面，计数
            if (line_cross_surface(prev, curr, thres, sam_num) < 0.999)
                count_surface_line++;
            if (i == psize - 2)
            {
                thres = min(all_nodes[path[i]].center_value, all_nodes[path[i + 1]].center_value);
                if (line_cross_surface(curr, next, thres, sam_num) < 0.999)
                    count_surface_line++;
            }

            double angle_deg = abs_angle(prev - curr, next - curr) / 180.0;
            angle_product *= angle_deg;
            angle_degrees.push_back(angle_deg);
        }
        double e_index = psize - 2 + 1.5 * count_inner - 1.0 * count_surface_line / (psize - 1);
        //double e_index = psize - 2 + 3*count_inner - 1.5* count_surface_line / (psize - 1);
        e_index = max(min(e_index, 9.0), 1.0);
        //if (e_index > 9) e_index = 9;
        //int e_index = psize - 2 + count_inner - count_surface_line;
        if (angle_product < 1e-8) translucency_score = 0.0;
        else  translucency_score = std::pow(angle_product, 1.0 / e_index);// / Kernels.size();

        if (show_debug)
        {
            cout << "angle_deg: ";
            for (auto ang_ : angle_degrees)
                cout << ang_ << "  ";
            cout << "  angle_deg product: " << angle_product << endl;
            cout << "count_inner: " << count_inner << "   count_surface_line: " << count_surface_line << "  e index:" << e_index << "    " << 1.0 / e_index << endl
                << translucency_score << "  * " << log(psize) << " = " << translucency_score * log(psize) << endl;
        }


        //cout << "Translucency_score2: " << translucency_score << "     "<< translucency_score * log(psize)<<endl;
        //translucency_score = std::pow(translucency_score, 1.0 / (psize - 2)) * (1+ 1.0* count_inner / psize);// / Kernels.size();
    }

    return translucency_score * log(psize);

}

double ModelGenerator::cal_kernel_translucency(int p_index, int & max_s1, int & max_s2, std::vector<int>& max_path, AdjacencyList adj, bool debug)  //计算点p所在的所有路径中通透性最大的一条路径作为通透性值，单独处理内部点且只有一条边以加速
{
    double ave_perm = 0.0;
    int count_ = 0;
    double max_perm = 0.0;
    max_s1 = -1;
    max_s2 = -1;
    max_path.clear();
    int kernel_num = Kernels.size();
    PathQuery p_bfs(kernel_num, adj, p_index);
    if (!p_bfs.isConnected)  //不连通，直接跳过
    {
        //cout << "有一条为空说明不再连通" << endl;
        return 0.0;
    }
	//if (adj[p_index].size() < 2 && (!Kernels[p_index].on_surface))    //内部点且只有一条连接边，无法形成路径，直接返回0
 //   {
	//	return max_perm;
 //   }
    // 双重循环遍历所有 s1, s2 组合, 复杂度 O(K^2)，其中 K 是 surface_points 的数量
    for (int i = 0; i < surface_kernels.size(); i++) 
    {
        // 剪枝：如果 s1 无法到达 p_index，则跳过
        int s1 = surface_kernels[i];
        for (int j = i+1; j < surface_kernels.size(); j++) 
        {
            int s2 = surface_kernels[j];
            double euclidean_dist = distance(Kernels[s1].center, Kernels[s2].center);
            if (debug)
                cout << "s1: " << s1 << "  s2: " << s2 << endl;
            std::vector<int> path1 = p_bfs.query_path(s1);
            std::vector<int> path2 = p_bfs.query_path(s2);

            std::vector<int> path_(path1.rbegin(), path1.rend());
            path_.insert(path_.end(), path2.begin() + 1, path2.end()); // 从 path2 的倒数第二个元素开始添加
            double graph_dist = cal_path_graph_length(path_);
			double path_translucency = calculate_path_translucency(path_);
            if(debug)
                cout << "s1: " << s1 << "  s2: " << s2 << endl
                << "angle trans: " << path_translucency << "   length ratio: " << euclidean_dist / graph_dist << endl;
            // 3. 计算通透性
            double perm = path_translucency;// *euclidean_dist / graph_dist;
            //ave_perm += perm;
            count_++;
            // 4. 更新最大值
            if (perm > max_perm)
            {
                max_perm = perm;
                max_s1 = s1;
                max_s2 = s2;
                max_path = path_;
            }
        }
    }

    return max_perm;
}

double ModelGenerator::cal_total_translucency(std::vector<GaussianKernel> gau,  AdjacencyList adj)
{
    int kernels_num = Kernels.size();
    double total_score = 0.0;
    max_paths_kernel.clear();
    kernel_translucency.clear();
    Paths.clear();
    //for (int i = 3; i < 4; i++)
    //for (int i = 53; i < 54; i++)
    int good_count = 0;
    for (int i = 0; i < kernels_num; i++)
    {
        int start = -1, end = -1;
        std::vector<int> max_path;
        double score_p = cal_kernel_translucency(i, start, end, max_path, adj, false);
        total_score += score_p;
        max_paths_kernel.push_back(make_pair(start, end));
		kernel_translucency.push_back(score_p);
        Paths.push_back(max_path);
        if (score_p >= Trans_thres)
            good_count++;
        if(standard_show)
        {
            cout << "Kernel " << i << " :  max translucency: " << score_p << "   from " << start << " to " << end << endl << "Max path: ";
			if (max_path.empty()) cout << " Empty path! ";
            else
                for (auto p : max_path) cout << p << "  ";
            cout << endl;
        }
    }
    cout << good_count << " Kernels have met the translucency threshold!" << endl;
    return total_score/ kernels_num;
}

vector<int> ModelGenerator::check_inner_leafs(vector<int> leafs_index)
{
    vector<int> inner_leafs;
    for (auto t : leafs_index)
        if (!Kernels[t].on_surface)
        {
            cout << "Inner Leaf Kernel: " << t << endl;
            inner_leafs.push_back(t);
        }
    int inner_leafs_num = inner_leafs.size();
    cout << "Inner_leafs_num: " << inner_leafs_num << endl;
    return inner_leafs;
}


int ModelGenerator::find_nearest_grid(Eigen::Vector3d point)
{
    int res = Resolution;
    double dx = 1.0 / (res - 1);
    // std::round 自动寻找最近的整数，即寻找最近的网格点
    int ix = static_cast<int>(std::round((point.x() + 0.5) / dx));
    int iy = static_cast<int>(std::round((point.y() + 0.5) / dx));
    int iz = static_cast<int>(std::round((point.z() + 0.5) / dx));
    //cout << "point: "<<point <<"  "<<ix<<"   "<<iy<<"  "<<iz<< endl;
    // 边界检查, 防止采样点略微超出 (0,0,0)~(1,1,1) 导致索引越界
    ix = std::max(0, std::min(ix, res - 1));
    iy = std::max(0, std::min(iy, res - 1));
    iz = std::max(0, std::min(iz, res - 1));

    int index = ix + iy * res + iz * res * res;

    return index;

}

double ModelGenerator::line_cross_surface(Eigen::Vector3d p1, Eigen::Vector3d p2, double thres, int sam_num)
{
    double sum_sdf = 0;
    for (int t = 1; t <= sam_num; t++)
    {
        double dt = 1.0 / (sam_num + 1);
        Eigen::Vector3d point = p1 * (1 - t * dt) + p2 * (t * dt);
        double sdf_val = SDF_ini[find_nearest_grid(point)];
        sum_sdf += sdf_val;
    }
    double center_sdf = sum_sdf / sam_num;
    if (center_sdf < -0.05)
    //if (center_sdf +1e-2 < thres )
    {
        //cout << "the line cross the model!" << endl;
        return 1.0; // 横贯模型
    }
    else
    {
        return max(1.0 - abs(center_sdf + 0.05) / 0.1, 0.0); //表面路径，无内部节点，通透性为1-dis
        //return 1.0 - abs((center_sdf - thres) / thres); //表面路径，无内部节点，通透性为1-dis
    }
}

double ModelGenerator::calculate_score(std::vector<std::vector<int>>  Paths)
{
    double sum_score = 0.0;
    for (auto t : Tube_edges)
    {
        sum_score += t.weight;
    }
	cout << "Total score: " << sum_score << endl;
	return sum_score;

    double total_weighted_permeability = 0.0;
    double total_path_length_sum = 0.0;
    if (Paths.empty()) {
        return 0.0;
    }
    else
    {
        for (auto path_ : Paths)
        {
            //std::pair<double, double>  properties = calculate_path_translucency(path_);
            //double path_permeability = properties.first;
            //double path_length = properties.second;

            //total_weighted_permeability += path_permeability * path_length;
            //total_path_length_sum += path_length;
        }
    }
    // --- 1. 计算通透性 Permeability(G) ---

    double permeability_G = 0.0;
    if (total_path_length_sum > 1e-9) {
        permeability_G = total_weighted_permeability / total_path_length_sum;
    }
    else {
        // 如果没有边界-边界路径（例如少于2个边界点），则通透性为0
        permeability_G = 0.0;
    }

    // --- 2. 计算成本 Cost(G) ---
    double cost_L = 0.0; // 长度成本
    int cost_E = Tube_edges.size();      // 边数成本

    for (auto& edge : Tube_edges)
    {
        cost_L += (Kernels[edge.from].center - Kernels[edge.to].center).norm();  //length
        //cost_L += edge.weight;  //length         
    }
     
	double w_e = 0.01; // 边数权重
    double cost_G = cost_L + w_e * cost_E;

    // --- 3. 计算最终得分 Score(G) ---
    if (cost_G < 1e-9) {
        // 成本几乎为0，如果通透性也为0，得分为0；否则得分非常高
        return (permeability_G > 1e-9) ? std::numeric_limits<double>::max() : 0.0;
    }

    return permeability_G / cost_G;
}


vector<int> ModelGenerator::cal_edge_usage(std::vector<std::vector<int>> Paths, bool show_debug)
{
    std::map<std::pair<int, int>, int> edge_count;
	vector<int> edge_usage_count;
    int count = 0;
    for(auto path : Paths)
    {
        //cout << "path " << count++ << "    : ";
        //for (auto p_in : path) cout << p_in << "  ";
        //cout << endl;
        if (path.size() < 2)
        {
			cout << "This path has less than 2 nodes." << endl;
            continue;
        }

        for (size_t i = 0; i < path.size() - 1; ++i) {
            int u = path[i];
            int v = path[i + 1];

            auto edge = std::minmax(u, v);
            edge_count[edge]++;
        }
    }
    //cout << "edge_count size: " << edge_count.size() << endl;
    for (int ei = 0; ei < Tube_edges.size(); ei++) 
    {
        Edge e = Tube_edges[ei];
        auto edge = std::minmax(e.from, e.to);
        edge_usage_count.push_back(edge_count[edge]);
        if(show_debug)
            cout << "Edge: " << ei <<"  from "<< e.from<<"  to "<<e.to <<"  is used "<< edge_count[edge] <<"  times"<<endl;
    }
    
    return edge_usage_count;
}

pair<double, double> ModelGenerator::add_edges(Edge cand_edge, AdjacencyList adj, std::vector<int>& max_path1, std::vector<int>& max_path2, bool debug)
{
    int kernels_num = Kernels.size();
    AdjacencyList adj_new = adj;
    //update adj_list
    if (cand_edge.from < kernels_num && cand_edge.to < kernels_num) {
        adj_new[cand_edge.from].push_back(cand_edge.to);
        adj_new[cand_edge.to].push_back(cand_edge.from);
    }

	int p1 = cand_edge.from;
	int p2 = cand_edge.to;
	int start = -1, end = -1;
    double length = length_path(p1, p2);
    double p1_add = cal_kernel_translucency(p1, start, end, max_path1, adj_new, false);
    p1_add = p1_add - kernel_translucency[p1];
    //show_path(max_path1);
    start = -1;
    end = -1;
    double p2_add = cal_kernel_translucency(p2, start, end, max_path2, adj_new);
	p2_add = p2_add - kernel_translucency[p2];

    /*double thres = min(Kernels[start_].center_value, Kernels[end_].center_value);
    translucency_score = line_cross_surface(Kernels[start_].center, Kernels[end_].center, thres, sam_num);*/

	//cout << "cal kernel score: " << p1_add << "  " << p2_add << "   cand_edge.length: "<< cand_edge.length<< "   delta_s: "<< (p1_add + p2_add) / cand_edge.length<<endl;
    //double score_add = (p1_add + p2_add)/ cand_edge.length;  //单位路径增加的通透性

    return make_pair(p1_add, p2_add);

}

bool ModelGenerator::replace_edges(int p_index, int replace_e, std::vector<Edge>& tube_edges, AdjacencyList& adj, AdjacencyList& unused_adj)
{
    int kernels_num = Kernels.size();
    std::vector<Edge> edge_mst_new = tube_edges;
    AdjacencyList adj_new = adj;
    AdjacencyList unused_adj_new = unused_adj;
	Edge re_edge = tube_edges[replace_e];
    //cout << "adj_new[" << p_index << "] size: " << adj_new[p_index].size() << "   ";
    //for (auto d : adj_new[p_index]) cout <<d << "   ";
    
    //update adj_list
    if (re_edge.from < kernels_num && re_edge.to < kernels_num) {
        auto& adj_a = adj_new[re_edge.from];
        adj_a.erase(std::remove(adj_a.begin(), adj_a.end(), re_edge.to), adj_a.end());

        auto& adj_b = adj_new[re_edge.to];
        adj_b.erase(std::remove(adj_b.begin(), adj_b.end(), re_edge.from), adj_b.end());
    }
    //这里暂时不修改unused_adj_new，不需要再计算这条边

    double max_delta_score = 0;
    int chosen_cand = -1;
    Edge chosen_e;
    pair<double, double> delta_score_max_pair;
    std::vector<int> max_path1, max_path2;
    //cout << "unused_adj_new[" << p_index << "] size: " << unused_adj_new[p_index].size() << "   ";
    //for (auto d : unused_adj_new[p_index]) cout << d << "   ";

    for (int candidate_p : unused_adj_new[p_index])
    {
        double dist = distance(Kernels[p_index].center, Kernels[candidate_p].center);
        double dist_w = dist * calculate_edge_weight(Kernels[p_index], Kernels[candidate_p]);
        Edge cand_edge = { p_index, candidate_p, dist, dist_w };
        std::vector<int> path1, path2;
        pair<double, double> delta_score = add_edges(cand_edge, adj_new, path1, path2, true);
        double score_add = (delta_score.first + delta_score.second) / cand_edge.length;  //单位路径增加的通透性
        if (delta_score.first > max_delta_score && delta_score.second > -1e-9)  //记录p_index通透值增大最多的边，且整体通透性要增加
        {
            max_delta_score = delta_score.first;
            chosen_cand = candidate_p;
            chosen_e = cand_edge;
            delta_score_max_pair = delta_score;
            max_path1 = path1;
            max_path2 = path2;
        }
    }
    if (chosen_cand != -1)
    {
        edge_mst_new[replace_e] = chosen_e; //替换
        //update adj_list
        if (p_index < Kernels.size() && chosen_cand < Kernels.size()) {
            adj_new[p_index].push_back(chosen_cand);
            adj_new[chosen_cand].push_back(p_index);
        }
        unused_adj_new[p_index].erase(std::remove(unused_adj_new[p_index].begin(), unused_adj_new[p_index].end(), chosen_cand), unused_adj_new[p_index].end());
        unused_adj_new[chosen_cand].erase(std::remove(unused_adj_new[chosen_cand].begin(), unused_adj_new[chosen_cand].end(), p_index), unused_adj_new[chosen_cand].end());
        unused_adj_new[re_edge.from].push_back(re_edge.to);  //删掉的边
        unused_adj_new[re_edge.to].push_back(re_edge.from);
        kernel_translucency[p_index] += delta_score_max_pair.first;
        kernel_translucency[chosen_cand] += delta_score_max_pair.second;
        Paths[p_index] = max_path1;
        Paths[chosen_cand] = max_path2;
 
        std::cout <<"Replace Edge ("<< re_edge.from<<", "<< re_edge.to<<") to Edge ("<< p_index << ", " << chosen_cand << "),  increasing Kernel " << p_index << " 's score to " << kernel_translucency[p_index] << "  by: " << delta_score_max_pair.first <<"and "<< delta_score_max_pair.second<< endl;
        tube_edges = edge_mst_new;
        adj = adj_new;
        unused_adj = unused_adj_new;
        return true;
    }
    else
    {
        std::cout << "No useful candidate edge found to replace!" << endl;
        return false;
    }

}


void ModelGenerator::optimize_mst(int opt_times_once, int edge_max, bool debug)
{
    double thres_tran = Trans_thres;
    int kernels_num = Kernels.size();
	int curr_edge_num = Tube_edges.size();
    int edge_add = edge_max * 0.1;
    int edge_max_num = max(edge_max, curr_edge_num+ edge_add);
    std::cout << "curr_edge: " << curr_edge_num << "   max_edge: " << edge_max_num << endl;

    int opt_count_total = 0;
    int num_add = 0, num_replace = 0;
	//check each kernel's translucency
    std::cout << "4.1 Adding edges for inner kernel with single adj:" << endl;
    for (int i = 0; i < kernels_num; i++)
    {
        //std::cout << "Adj_list[i].size() " << Adj_list[i].size() << endl;
        if (Adj_list[i].size() < 2 && (!Kernels[i].on_surface)) //内部点且只有一条连接边，需要直接加边
        {
            std::cout << "===========================================================" << endl<<
                "Kernel " << i << " has a low translucency: " << kernel_translucency[i] << "  lower than thres_tran:" << thres_tran << endl;
            double max_delta_score = 0;
            int chosen_cand = -1;
            Edge chosen_e;
            pair<double, double> delta_score_max_pair;
            std::vector<int> max_path1, max_path2;
            for (int candidate_p : Unused_adj_list[i])
            {
                double dist = distance(Kernels[i].center, Kernels[candidate_p].center);
                double dist_w = dist * calculate_edge_weight(Kernels[i], Kernels[candidate_p]);
                Edge cand_edge = { i, candidate_p, dist, dist_w };
                std::vector<int> path1, path2;
                pair<double, double> delta_score = add_edges(cand_edge, Adj_list, path1, path2);
                double score_add = (delta_score.first + delta_score.second) / cand_edge.length;  //单位路径增加的通透性
                if (delta_score.first > max_delta_score && delta_score.second > -1e-9)
                {
                    max_delta_score = delta_score.first;
                    chosen_cand = candidate_p;
                    chosen_e = cand_edge;
                    delta_score_max_pair = delta_score;
                    max_path1 = path1;
                    max_path2 = path2;
                }
            }
            if (chosen_cand != -1)
            {
                Tube_edges.push_back(chosen_e);
                //update adj_list
                if (i < Kernels.size() && chosen_cand < Kernels.size()) {
                    Adj_list[i].push_back(chosen_cand);
                    Adj_list[chosen_cand].push_back(i);
                }
                Unused_adj_list[i].erase(std::remove(Unused_adj_list[i].begin(), Unused_adj_list[i].end(), chosen_cand), Unused_adj_list[i].end());
                Unused_adj_list[chosen_cand].erase(std::remove(Unused_adj_list[chosen_cand].begin(), Unused_adj_list[chosen_cand].end(), i), Unused_adj_list[chosen_cand].end());
                kernel_translucency[i] += delta_score_max_pair.first;
                kernel_translucency[chosen_cand] += delta_score_max_pair.second;
                Paths[i] = max_path1;
                Paths[chosen_cand] = max_path2;
                curr_edge_num++;
                opt_count_total++;
                num_add++;
                cout << "Add edge from " << i << " to " << chosen_cand << "  increasing Kernel " << i << " 's score to " << kernel_translucency[i] << "  by: " << delta_score_max_pair.first << endl << "Edge length: " << chosen_e.length << "   new edge num : " << curr_edge_num <<  endl << endl;
            }
        }
    }

    cout << num_add << " edges have beed added in this step!" << endl;

    double finalTranslucency = cal_total_translucency(Kernels, Adj_list);
    std::cout << "After optimization, total score increases to " << finalTranslucency << " with edges to " << Tube_edges.size() << endl;
    cal_edge_usage(Paths, NO_DEBUG);

    cout << endl<<"4.2 Iteratively optimizing edges:" << endl;
    num_add = 0;
    for(int i = 0; i < kernels_num; i++)
    {
        int opt_count = 0;
        std::cout << "Kernel " << i << " translucency: " << kernel_translucency[i] << std::endl;
        if (kernel_translucency[i] < 1e-9 && (!Kernels[i].on_surface)) //内部点且没有通路，有可能子节点已经连通，所以首先重新计算通透性再判断是否修改
        {
			int start = -1, end = -1;
            std::vector<int> max_path;
            kernel_translucency[i] = cal_kernel_translucency(i, start, end, max_path, Adj_list, false);
            
            if (kernel_translucency[i] > 1e-9)
                cout << "Kernel: "<<i <<"'s translucency is increased to " << kernel_translucency[i] << " due to the last optimization!" << endl;
            if (kernel_translucency[i] >= thres_tran)
            {
                Paths[i] = max_path;
            }
        }
		//接下来是循环优化
        bool add_end = (curr_edge_num < edge_max_num) ? false : true;
        bool replace_end = false;
        while(kernel_translucency[i] < thres_tran && opt_count < opt_times_once && opt_count_total< opt_times_once * kernels_num)
        {
            opt_count++;
            //std::cout << "===========================================================" << endl <<
            //    "Kernel " << i << " with translucency: " << kernel_translucency[i] << " < " << thres_tran << endl;

            if (add_end && replace_end)
            {
                //cout << "Both Adding and Replacing END for Kernel " << i << "! Exit optimization for this kernel!" << endl;
				break;
            }
            if (curr_edge_num < edge_max_num && !add_end)  //没有到边数上限且新增边能增加score，则选择添加边
            {
                double max_delta_score = 0;
                int chosen_cand = -1;
                Edge chosen_e;
                pair<double, double> delta_score_max_pair;
                std::vector<int> max_path1, max_path2;

                for (int candidate_p : Unused_adj_list[i])
                {
                    double dist = distance(Kernels[i].center, Kernels[candidate_p].center);
                    double dist_w = dist * calculate_edge_weight(Kernels[i], Kernels[candidate_p]);
                    Edge cand_edge = { i, candidate_p, dist, dist_w };
                    std::vector<int> path1, path2;
                    pair<double, double> delta_score = add_edges(cand_edge, Adj_list, path1, path2);
                    double score_add = (delta_score.first + delta_score.second) / cand_edge.length;  //单位路径增加的通透性
                  
                    if (delta_score.first > max_delta_score && delta_score.second > -1e-9)
                    {
                        max_delta_score = delta_score.first;
                        chosen_cand = candidate_p;
                        chosen_e = cand_edge;
                        delta_score_max_pair = delta_score;
                        max_path1 = path1;
						max_path2 = path2;
                    }
                }
                if (chosen_cand != -1)
                {
                    Tube_edges.push_back(chosen_e);
                    //update adj_list
                    if (i < Kernels.size() && chosen_cand < Kernels.size()) {
                        Adj_list[i].push_back(chosen_cand);
                        Adj_list[chosen_cand].push_back(i);
                    }
					Unused_adj_list[i].erase(std::remove(Unused_adj_list[i].begin(), Unused_adj_list[i].end(), chosen_cand), Unused_adj_list[i].end());
                    Unused_adj_list[chosen_cand].erase(std::remove(Unused_adj_list[chosen_cand].begin(), Unused_adj_list[chosen_cand].end(), i), Unused_adj_list[chosen_cand].end());
                    kernel_translucency[i] += delta_score_max_pair.first;
                    kernel_translucency[chosen_cand] += delta_score_max_pair.second;
					Paths[i] = max_path1;
					Paths[chosen_cand] = max_path2;
                    curr_edge_num++;
                    opt_count_total++;
                    num_add++;
                    //cout <<"Add edge from " << i << " to " << chosen_cand << "  increasing Kernel "<<i<<" 's score to "<< kernel_translucency[i] <<"  by: "<< delta_score_max_pair.first <<endl<<"Edge length: "<< chosen_e.length<<"   new edge num : " << curr_edge_num << endl << endl;
                }
                else
                {
                    //cout << "No useful candidate edge to add!" << endl;
                    add_end = true;
                }
            }
            else if(!replace_end) //到达边数上限，或者增加边无法提升得分，则选择替换边
            {
                //if(curr_edge_num >= edge_max_num) 
                //    cout << "Reach the edge LIMITATION, REPLACE the edge!" << endl;
                //else if(add_end)
                //    cout<< "Adding END, then REPLACE the edge!" << endl;

                vector<int> edge_importance = cal_edge_usage(Paths, NO_DEBUG);
                std::vector<std::pair<int, int>> sorted_edges;
                sorted_edges.reserve(edge_importance.size());
                for (int ei = 0; ei < edge_importance.size(); ++ei) {
                    sorted_edges.push_back({ ei, edge_importance[ei]});
                }
                std::sort(sorted_edges.begin(), sorted_edges.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                    // 如果次数不同，按次数从小到大排序 (降序)
                    if (a.second != b.second) {
                        return a.second < b.second;
                    }
                    return a.first < b.first; // 如果次数相同，按索引从小到大排序 (保持稳定性，可选)
                    });
                //cout << "sorted_edges.size(): " << sorted_edges.size() << endl;
				bool rep_get = false;
                for (int s_index = 0; s_index < sorted_edges.size(); s_index++) 
                {
                    std::pair<int, int> sorted_e = sorted_edges[s_index];
					int replace_e = sorted_e.first;
                    if ( i != Tube_edges[replace_e].from && i != Tube_edges[replace_e].to)  //非kernel i 的边
                    {
                        //cout << "The edge ("<< Tube_edges[replace_e].from<<", "<< Tube_edges[replace_e].to<<") is not one of Kernel i!Skip!" << endl;
                        continue;
                    }
                    else  if (sorted_e.second != 0 && !find_edge_in_path(Tube_edges[replace_e], Paths[i])) //有用过，但不在kernel i 的最大通透性路径上，跳过
                    {
                        //cout << "The edge (" << Tube_edges[replace_e].from << ", " << Tube_edges[replace_e].to << ") is used in other paths!  Skip!" << endl;
                        continue;
                    }
                    else if (sorted_e.second > 1) //使用次数较多，跳过
                    {
                        //cout << "The edge (" << Tube_edges[replace_e].from << ", " << Tube_edges[replace_e].to << ") is used " << sorted_e.second << " times!  Skip!" << endl;
                        continue;
					}
					else  
                    {
                        //尝试替换这条边
                        //Edge removed_edge = Tube_edges[replace_i];
                        //cout << "Optimization Kernel "<<i << ": Edge "<< replace_e <<"("<< Tube_edges[replace_e].from<<", "<< Tube_edges[replace_e].to<<") with usage count : " << sorted_e.second <<" will be replaced ... " << endl;
                        if (Adj_list[i].size() < 2) //虽然替换，但是只有一条连接边，无法替换
                        {
                            //cout << "improve the edge limitation" << endl;
                            edge_max_num++;
                            break; //提升上限，跳出循环
                        }
                        else if (replace_edges(i, replace_e, Tube_edges, Adj_list, Unused_adj_list))
                        {
                            opt_count_total++;
							num_replace++;
                            //cout << "The num of " << num_replace << " edges have beed repalced!" << endl << endl;
                            rep_get = true;
                            break; //成功替换，跳出循环
						}
                    }
                }
                if (!rep_get) //说明替换遍历完了
                    replace_end = true;
            }
        }
        if(opt_count >=  opt_times_once)
            //cout << "-------------------------------- Reach one edge's limitation! -------------------------------- " << endl
            //<<"In optimization iteration of kernel " <<i<<": "<< num_add << " edges added and "<<num_replace<<" edges replaced!" << endl;;
        if (opt_count_total >= opt_times_once * kernels_num)
        {
            //cout << "--------------------------------  Reach the total optimization limitation! -------------------------------- " << endl;
            break;
        }
            
    }
    cout << "During the whole optimization, " << opt_count_total << " edge, " << "including " << num_replace << " repalced and " << opt_count_total - num_replace << " added, have beed added and repalced!" << endl << endl;
}



void ModelGenerator::optimize_mst2(int itea_max_times, int max_edge, bool iter_add, bool debug) //max_edge = 0代表最大边数递增，!=0代表固定最大边数
{
    std::cout << "--------------------4. Optimizing the connection trees --------------------" << endl;
	int origin_edge_num = Tube_edges.size();
    int curr_edge_num = origin_edge_num;
	int ratio = 1;
	double delta_trans_score = 1000.0;
    double thres_tran = Trans_thres;
    int kernels_num = Kernels.size();
    int edge_max_num;
    int opt_count = 0, opt_add = 0, opt_replace = 0;
    double new_trans_score = finalTranslucency;

    while (ratio < itea_max_times && delta_trans_score>0.01)
    {
        if(iter_add)
            edge_max_num = min(max_edge, (int)(origin_edge_num * (1.1 + ratio / 10.0)));
		else
            edge_max_num = max_edge;

        vector<pair<int, double>> kernel_info_order;
        for (int k = 0; k < kernel_translucency.size(); k++)
            kernel_info_order.push_back(make_pair(k, kernel_translucency[k])); 
        sort_min2max(kernel_info_order);
        //for (int i = 0; i < kernels_num; i++)
        for (int index = 0; index < kernels_num; index++)
        {
            int i = kernel_info_order[index].first;
            double origin_trans = kernel_translucency[i];
            int start = -1, end = -1;
            std::vector<int> max_path;
            kernel_translucency[i] = cal_kernel_translucency(i, start, end, max_path, Adj_list, debug);
            if (origin_trans != kernel_translucency[i])
            {
                Paths[i] = max_path;
                cout << "Kernel " << i << " translucency has been changed to: " << kernel_translucency[i]<<" by "<< kernel_translucency[i] - origin_trans << " due to the last optimization!" << endl;
            }
            else
                std::cout << "Kernel " << i << " translucency: " << kernel_translucency[i] << std::endl;

            if (kernel_translucency[i] > thres_tran)  //符合要求
            {
                cout << "Good Kernel"<<endl;
                continue;
            }
            
            double delta_add_score = 0.0, delta_replace_score = 0.0;
            bool choose_add = false;
			bool choose_replace = false;
            int chosen_cand = -1;
            int chosen_cand_r = -1;
            int replaced_edge_index = -1;
            Edge chosen_e, chosen_e2;
            pair<double, double> delta_score_max_pair, delta_score_max_pair2;
            std::vector<int> max_path1, max_path2, max_path3, max_path4;

            //-------------------calculate add edge score increase---------------------------
            for (int candidate_p : Unused_adj_list[i])
            {
                double dist = distance(Kernels[i].center, Kernels[candidate_p].center);
                double dist_w = dist * calculate_edge_weight(Kernels[i], Kernels[candidate_p]);
                Edge cand_edge = { i, candidate_p, dist, dist_w };
                std::vector<int> path1, path2;
                pair<double, double> delta_score = add_edges(cand_edge, Adj_list, path1, path2);
                if (delta_score.first > delta_add_score && delta_score.second > -1e-9)
                {
                    delta_add_score = delta_score.first;
                    chosen_cand = candidate_p;
                    chosen_e = cand_edge;
                    delta_score_max_pair = delta_score;
                    max_path1 = path1;
                    max_path2 = path2;
                }
            }

            if (Adj_list[i].size() >= 2 || Kernels[i].on_surface) //外部点或者有多条连接边，计算替换
            {
                //-------------------calculate add edge score increase---------------------------
                vector<int> edge_importance = cal_edge_usage(Paths, NO_DEBUG);
                std::vector<std::pair<int, int>> sorted_edges;
                for (int ei = 0; ei < edge_importance.size(); ++ei) {
                    sorted_edges.push_back({ ei, edge_importance[ei] });
                }
                sort_min2max(sorted_edges);// 如果次数不同，按次数从小到大排序 (降序), 如果 second 相等，希望保持原顺序

                for (int s_index = 0; s_index < sorted_edges.size(); s_index++)
                {
                    std::pair<int, int> sorted_e = sorted_edges[s_index];
                    int replace_e = sorted_e.first;
                    if (i != Tube_edges[replace_e].from && i != Tube_edges[replace_e].to) {  //非kernel i 的边
                        //cout << "The edge ("<< Tube_edges[replace_e].from<<", "<< Tube_edges[replace_e].to<<") is not one of Kernel i!Skip!" << endl;
                        continue;
                    }
                    else  if (sorted_e.second != 0 && !find_edge_in_path(Tube_edges[replace_e], Paths[i])) { //有用过，但不在kernel i 的最大通透性路径上，跳过
                        //cout << "The edge (" << Tube_edges[replace_e].from << ", " << Tube_edges[replace_e].to << ") is used in other paths!  Skip!" << endl;
                        continue;
                    }
                    else if (sorted_e.second > 1) { //使用次数较多，跳过
                        //cout << "The edge (" << Tube_edges[replace_e].from << ", " << Tube_edges[replace_e].to << ") is used " << sorted_e.second << " times!  Skip!" << endl;
                        continue;
                    }
                    else
                    {
                        //尝试替换这条边
                        std::cout << "Optimization Kernel " << i << ": Edge " << replace_e << "(" << Tube_edges[replace_e].from << ", " << Tube_edges[replace_e].to << ") with usage count : " << sorted_e.second << " will be replaced ... " << endl;
                        std::vector<Edge> edge_mst_new = Tube_edges;
                        AdjacencyList adj_new = Adj_list;
                        AdjacencyList unused_adj_new = Unused_adj_list;
                        Edge re_edge = edge_mst_new[replace_e];

                        if (re_edge.from < kernels_num && re_edge.to < kernels_num) {
                            auto& adj_a = adj_new[re_edge.from];
                            adj_a.erase(std::remove(adj_a.begin(), adj_a.end(), re_edge.to), adj_a.end());

                            auto& adj_b = adj_new[re_edge.to];
                            adj_b.erase(std::remove(adj_b.begin(), adj_b.end(), re_edge.from), adj_b.end());
                        }
                        //这里暂时不修改unused_adj_new，不需要再计算这条边

                        for (int candidate_p : unused_adj_new[i])
                        {
                            double dist = distance(Kernels[i].center, Kernels[candidate_p].center);
                            double dist_w = dist * calculate_edge_weight(Kernels[i], Kernels[candidate_p]);
                            Edge cand_edge = { i, candidate_p, dist, dist_w };
                            std::vector<int> path1, path2;
                            pair<double, double> delta_score = add_edges(cand_edge, adj_new, path1, path2, debug);
                            if (delta_score.first > delta_replace_score && delta_score.second > -1e-9)  //记录p_index通透值增大最多的边，且整体通透性要增加
                            {
                                delta_replace_score = delta_score.first;
                                chosen_cand_r = candidate_p;
                                chosen_e2 = cand_edge;
                                delta_score_max_pair2 = delta_score;
                                max_path3 = path1;
                                max_path4 = path2;
                                replaced_edge_index = replace_e;
                            }
                        }
                    }
                }
                //compare delta_add_score and delta_replace_score
                if (chosen_cand != -1 && delta_add_score >= delta_replace_score && curr_edge_num <= edge_max_num)
                    choose_add = true;
                else if (chosen_cand_r != -1 && delta_add_score < delta_replace_score)
                    choose_replace = true;
                else
                    std::cout << "No LEGAL operation found for kernel " << i << "!" << std::endl;
                cout<< "chosen_cand: " << chosen_cand << "  " << chosen_cand_r
                    << "     delta_score: " << delta_add_score << "   " << delta_replace_score
                    << "     edges: " << curr_edge_num << "  " << edge_max_num << endl << endl;

            }
            else  //内部点且只有一条连接边，需要直接加边
            {
                if (chosen_cand != -1)
                {
                    std::cout << "This inner Kernel with one adj needs another edge!" << endl;
                    choose_add = true;
                }
            }
			//execute the chosen operation
            if (choose_add)
            {
                Tube_edges.push_back(chosen_e);
                //update adj_list
                if (i < Kernels.size() && chosen_cand < Kernels.size()) {
                    Adj_list[i].push_back(chosen_cand);
                    Adj_list[chosen_cand].push_back(i);
                }
                Unused_adj_list[i].erase(std::remove(Unused_adj_list[i].begin(), Unused_adj_list[i].end(), chosen_cand), Unused_adj_list[i].end());
                Unused_adj_list[chosen_cand].erase(std::remove(Unused_adj_list[chosen_cand].begin(), Unused_adj_list[chosen_cand].end(), i), Unused_adj_list[chosen_cand].end());
                kernel_translucency[i] += delta_score_max_pair.first;
                kernel_translucency[chosen_cand] += delta_score_max_pair.second;
                Paths[i] = max_path1;
                Paths[chosen_cand] = max_path2;
                curr_edge_num++;
                opt_count++;
                opt_add++;
                std::cout << "Add edge from " << i << " to " << chosen_cand << "  increasing Kernel " << i << " 's score to " << kernel_translucency[i] << "  by: " << delta_score_max_pair.first << "and "<< delta_score_max_pair.second<<endl <<"new edge num : " << curr_edge_num << endl << endl;
            }
            else if (choose_replace)
            {
                Edge re_edge = Tube_edges[replaced_edge_index];
                if (re_edge.from < kernels_num && re_edge.to < kernels_num) {
                    auto& adj_a = Adj_list[re_edge.from];
                    adj_a.erase(std::remove(adj_a.begin(), adj_a.end(), re_edge.to), adj_a.end());

                    auto& adj_b = Adj_list[re_edge.to];
                    adj_b.erase(std::remove(adj_b.begin(), adj_b.end(), re_edge.from), adj_b.end());
                }
                Tube_edges[replaced_edge_index] = chosen_e2; //替换
                //update adj_list
                if (i < Kernels.size() && chosen_cand_r < Kernels.size()) {
                    Adj_list[i].push_back(chosen_cand_r);
                    Adj_list[chosen_cand_r].push_back(i);
                }
                Unused_adj_list[i].erase(std::remove(Unused_adj_list[i].begin(), Unused_adj_list[i].end(), chosen_cand_r), Unused_adj_list[i].end());
                Unused_adj_list[chosen_cand_r].erase(std::remove(Unused_adj_list[chosen_cand_r].begin(), Unused_adj_list[chosen_cand_r].end(), i), Unused_adj_list[chosen_cand_r].end());
                Unused_adj_list[re_edge.from].push_back(re_edge.to);  //删掉的边
                Unused_adj_list[re_edge.to].push_back(re_edge.from);
                kernel_translucency[i] += delta_score_max_pair2.first;
                kernel_translucency[chosen_cand_r] += delta_score_max_pair2.second;
                Paths[i] = max_path3;
                Paths[chosen_cand_r] = max_path4;
                opt_count++;
                opt_replace++;
                std::cout << "Replace Edge (" << re_edge.from << ", " << re_edge.to << ") to Edge (" << i << ", " << chosen_cand_r << "),  increasing Kernel " << i << " 's score to " << kernel_translucency[i] << "  by: " << delta_score_max_pair2.first << "and " << delta_score_max_pair2.second << endl<<endl;

            }
        }
        double last_trans_score = new_trans_score;
        new_trans_score = cal_total_translucency(Kernels, Adj_list);
        std::cout << "======================================After optimization "<< ratio++<<", total score increases from "<< last_trans_score <<" to " << new_trans_score << " with edges to " << Tube_edges.size() << "======================================"<<endl;
        delta_trans_score = new_trans_score - last_trans_score;
    }
	cout << "Finally, total iterations: " << ratio - 1<<", last delta: "<< delta_trans_score<<", optimization: "<<opt_count << ", including " << opt_add << " added and " << opt_replace << " replaced!" << endl;
}

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

    // 如果投影点不在线段上，返回正值（不在管道内）
    if (t_proj <= 0.0 || t_proj >= 1.0) {
        //return abs(abs(t_proj - 0.5) - 0.5);
        return 1;
    }
    // 将t限制在[0, 1]区间，这样空间中所有的点都会被映射到线段上最近的点
    t_proj = std::max(0.0, std::min(1.0, t_proj));
    Eigen::Vector3d p_proj = c1 + t_proj * line_dir; // p在轴线上的投影点

    //double r0 = sqrt(mid_radius_factor) * k1.sigma;
    //double r4 = sqrt(mid_radius_factor) * k2.sigma;
    double r0 = sqrt(mid_radius_factor) * k1.sigma_perp;
    double r4 = sqrt(mid_radius_factor) * k2.sigma_perp;

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
    //return iso_level_C - energy;
    return iso_level_C * (1 - std::exp(1.0 - energy));

}


double ModelGenerator::generate_tube2(Eigen::Vector3d& p, GaussianKernel& k1, GaussianKernel& k2, double iso_level_C, double mid_radius_factor)
{
    const Eigen::Vector3d& c1 = k1.center;
    const Eigen::Vector3d& c2 = k2.center;

    // 1. 计算 ω_i 和 ω_j, 这里根据sigma进行转换
// 假设 ω = 1/(2 * sigma^2)，保持与标准高斯函数一致
    //double omega1 = 1.0 / (2.0 * k1.sigma * k1.sigma);
    //double omega2 = 1.0 / (2.0 * k2.sigma * k2.sigma);
    double omega1 = 1.0 / (2.0 * k1.sigma_parallel * k1.sigma_perp);
    double omega2 = 1.0 / (2.0 * k2.sigma_parallel * k2.sigma_perp);
    double avg_omega = (omega1 + omega2) / 2.0;
    double mu = 10 * mid_radius_factor * avg_omega;  //mid_radius_factor越大，通道越细

    // 2. 计算点到线段的距离（论文中的 ||p - s_ij||）
    Eigen::Vector3d line_dir = c2 - c1;
    double line_length = line_dir.squaredNorm();

    Eigen::Vector3d w = p - c1;
    double dis;
    double dis_c = w.dot(line_dir);
    if (dis_c <= 0) {
        dis = w.squaredNorm();;
    }
    else if (dis_c >= line_length) {
        dis = (p - c2).squaredNorm();
    }
    else {
        double t_proj = dis_c / line_length;
        Eigen::Vector3d p_proj = c1 + t_proj * line_dir; // p在轴线上的投影点
        dis = (p - p_proj).squaredNorm();
    }

    // 3. 计算管道函数的第一部分：沿线段的高斯函数
    double tunnelMain = std::exp(-mu * dis);

    // 4. 计算要减去的两个孔隙函数部分
    double pore1 = k1.gaussian_fun(p);
    double pore2 = k2.gaussian_fun(p);
    // 5. 组合管道函数
    double tubeValue = tunnelMain +pore1 + pore2;
    //cout << "tunnelMain:  " << tunnelMain << "   " << pore1 << "   " << pore2 << endl;
    //return  tubeValue;
    return iso_level_C - tubeValue;
}

int ModelGenerator::generate_mst_tubes(int grid_num, int res, double iso, double gaus_iso, double smooth_t)
{
    //单独保存高斯孔隙场SDF
    Eigen::VectorXd SDF_gaussian(grid_num);
    Eigen::VectorXd SDF_gaussian_tubes(grid_num);
    SDF_out.resize(grid_num);
    double void_count = 0;
    double tube_radius = Tube_radius_factor;

    for (int idx = 0; idx < grid_num; ++idx) {
        Eigen::Vector3d p = GV.row(idx);
        //gaussian kernel
        SDF_gaussian(idx) = combinedSDF(p, Kernels, gaus_iso);
        //tubes
        double sdf_p = 1000.0;
        for (auto& e : Tube_edges) {
            sdf_p = min(sdf_p, generate_tube2(p, Kernels[e.from], Kernels[e.to], gaus_iso, tube_radius));
            //sdf_p = min(sdf_p, generate_tube(p, Kernels[e.from], Kernels[e.to], gauss_iso, tube_radius));
        }
        SDF_gaussian_tubes(idx) = smooth_UnionSDF(SDF_gaussian(idx), sdf_p, 3 * smooth_t);
    }

    if (Enable_noise)
    {
        double noise_amplitude = 0.1;
        double band_width = 0.08;
        double spatial_frequency = 1.8;
        double noise_level_ratio = 4.5;
        add_noise_near_isosurface(SDF_ini, GV, Isolevel, noise_amplitude, band_width, spatial_frequency);
        //add_noise_near_isosurface(SDF_gaussian_tubes, GV, Isolevel, 0.66, 0.2, 6.3);  //200分辨率
        //add_noise_near_isosurface(SDF_gaussian_tubes, GV, Isolevel, 0.8, 0.3, 7);  //估计300分辨率
        add_noise_near_isosurface(SDF_gaussian_tubes, GV, Isolevel, 0.58, 0.16, 10);  //100分辨率
    }

    for (int idx = 0; idx < grid_num; ++idx) {
        Eigen::Vector3d p = GV.row(idx);
        SDF_out(idx) = smooth_IntersecSDF(SDF_ini(idx), -SDF_gaussian_tubes(idx), smooth_t);
        if (SDF_out(idx) < iso) {
            void_count += 1;
        }
    }
    //std::cout << "成功在仿生形状内放置 " << void_centers.size() << " 个空洞点" << std::endl;

    // Marching Cubes
    MarchingCubes(SDF_out, GV, res, res, res, iso, V_out, F_out);   //final result

    Eigen::MatrixXd V_g; //输出网格顶点
    Eigen::MatrixXi F_g; // 输出网格面片
    MarchingCubes(SDF_gaussian, GV, res, res, res, iso, V_g, F_g);  //gaussian combined
    view_model(V_g, F_g);

    MarchingCubes(SDF_gaussian_tubes, GV, res, res, res, iso, V_t, F_t);  //gaussian combined with tubes
    //view_model(V_t, F_t);

    //view_three_models(V_out, F_out, V_t, F_t, V_t, F_t, Eigen::RowVector3d(1, 0, 0));

    std::cout << "Generated mesh: " << V_out.rows() << " vertices, " << F_out.rows() << " faces" << std::endl;

    return void_count;
}

void ModelGenerator::test_item()
{
    //std::cout << "Model A loaded successfully." << std::endl;
    Eigen::MatrixXd V2;
    Eigen::MatrixXi F2;
    scale_factor = Mesh2SDF(V_ini, F_ini, GV, SDF_ini, bb_min, bb_max);
    int res = Resolution;
    std::string raw_file = "Tshape.voi";
    saveSDFtoVOI(raw_file, SDF_ini, res, res, res);

    VoxelGrid grids = SDFtoVoxel(SDF_ini, bb_min, bb_max, res, res, res);
    SupportCheckResult scr = checkSupportVoxel(grids, 0.5);
    
    /*std::string npy_filename = "D:/VSprojects/TaihuStone/model/npy/"+ input_file+ "_voxelized_model_" + std::to_string(res) + "x" + std::to_string(res) + "x" + std::to_string(res) + ".npy";
    saveVoxelGridAsNPY(grids.rho, res, npy_filename);*/
    
    Eigen::VectorXd scalar(res * res * res);

    for (int i = 0; i < res * res * res; ++i) {
        //scalar(i) = -static_cast<double>(voxel_grid[i]);
        scalar(i) = -grids.rho[i];
    }

    igl::marching_cubes(scalar, GV, res, res, res, -0.5, V2, F2);
    V2 = V2 / scale_factor;
    view_model(V2, F2);

    igl::write_triangle_mesh("result_change.stl", V2, F2);


}


