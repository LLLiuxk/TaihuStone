#include "modelGen.h"

GaussianKernel::GaussianKernel(Eigen::Vector3d center_, double sigma_, double amplitude_, double center_value_)
{
    center = center_;
    sigma = sigma_;         // 核的大小/影响力范围 (高斯函数的标准差)
    amplitude = amplitude_;
    center_value = center_value_;
    on_surface = is_on_surface();
}

double GaussianKernel::gaussian_fun(const Eigen::Vector3d& p)
{
	double d2 = (p - center).squaredNorm();
	double gau_value = amplitude * std::exp(-d2 / (2.0 * sigma * sigma));   //sigma越大，函数越平缓
    //std::cout << "  gau_value:"<< gau_value <<std::endl;
	return gau_value;

}

bool GaussianKernel::is_on_surface() const
{
    double ratio = Gauss_level / amplitude;
    double d_value = sqrt(-2.0 * sigma * sigma * std::log(ratio));
    //cout << d_value + 1e-3 << "   " << abs(center_value) << endl;
    return d_value + 1e-3 > abs(center_value);
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
    std::vector<double> pore_sdfs;
    std::vector<int> inside_indices;
    sample_interior_points(pore_centers, pore_sdfs, inside_indices, pores, gen);
    double model_count = inside_indices.size();
    //for (int i = 0; i < pore_centers.size(); i++)   std::cout << "i: " << i << "  " << pore_centers[i] << std::endl;

    generate_gaussians(pore_centers, pore_sdfs, gen);

	int kernel_num = Kernels.size();
    double void_count = 0;
    SDF_out.resize(grid_num);
	//单独保存高斯孔隙场SDF
    Eigen::VectorXd SDF_gaussian(grid_num);
    Eigen::VectorXd SDF_gaussian_tubes(grid_num);
    int degree_limit = (kernel_num - 1) / (kernel_num - surface_kernels.size());
    degree_limit = max(degree_limit, Min_degree);
    Tube_edges = pores_connection_mst(Kernels, degree_limit);

	// 构建邻接表
    Adj_list = construct_adj_list(Tube_edges, kernel_num);


    vector<int> leafs_index = all_leafs_mst(Tube_edges);
    vector<int> inner_leafs;
    for (auto t : leafs_index)
        if (!Kernels[t].on_surface)
        {
            cout << "inner t: " << t << endl;
            inner_leafs.push_back(t);
        }
        
	cout << "inner_leafs_num: " << inner_leafs.size() << endl;
    //-----------test---------------
    
    double trans_score = cal_total_translucency(Kernels, surface_kernels);

    cout << "total score: " << trans_score << endl;

 //   int p_index = 8;
	//int max_s1 = -1, max_s2 = -1;
 //  double p_score = cal_kernel_translucency(p_index, max_s1, max_s2);
	////cout << "p_index: " << p_index << "   max_s1: " << max_s1 << "   max_s2: " << max_s2 << "   p_score: " << p_score << endl;

 //   std::vector<int> path_ = find_specified_path(0, 0, 8);
 //   double path_translucency = calculate_path_translucency(path_);
    
	//cout << "path_score: " << path_score << endl;
// ----------test construct_dist_map--------------------------------
 //   std::vector<int>  path_ = find_path_in_tree(0, 13, kernel_num);
	//for (auto p : path_) cout << "path_p: " << p << " ";
 //   double lg1 = length_graph_path(0, 13);
 //   std::vector<double> dis_map = construct_dist_map(0, Adj_list);
 //   cout << "lg1 : " << lg1 << "   lg2 : " << dis_map[13] << endl;


    std::vector<std::vector<int>>  surface_paths;
    for(auto k1: surface_kernels)
        for (auto k2 : surface_kernels)
            if (k1 != k2)
            {
                std::vector<int>  path_ = find_path_in_tree(k1, k2, kernel_num);
				surface_paths.push_back(path_);
            }
	std::cout << "Surface kernel paths number: " << surface_paths.size() << std::endl;

    


    /*double score_ = calculate_score(surface_paths);
    cout << "Tube_edges.size: " << Tube_edges.size() << "  score: " << score_ << endl;*/
    //generate tubes
	double tube_radius = Tube_radius_factor;
    for (int idx = 0; idx < grid_num; ++idx) {
        Eigen::Vector3d p = GV.row(idx);
        //gaussian kernel
		SDF_gaussian(idx) = combinedSDF(p, Kernels, gauss_iso);
		//tubes
        double sdf_p = 1000.0;
        for (auto& e : Tube_edges){
            sdf_p = min(sdf_p, generate_tube2(p, Kernels[e.from], Kernels[e.to], gauss_iso, tube_radius));
        }
        //for (int kx = 0; kx < Kernels.size() - 1; kx++)
        //{
        //    //sdf_p = min(sdf_p, generate_tube(p, Kernels[kx], Kernels[kx + 1], gauss_combined, tube_radius));
        //    sdf_p = min(sdf_p, generate_tube2(p, Kernels[kx], Kernels[kx + 1], gauss_iso, tube_radius));
        //}
        
        SDF_gaussian_tubes(idx) = smooth_UnionSDF(SDF_gaussian(idx), sdf_p, smooth_t);

        SDF_out(idx) = smooth_IntersecSDF(SDF_ini(idx), -SDF_gaussian_tubes(idx), smooth_t);
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
    

    Eigen::VectorXd SDF_gaussian_tubes2(grid_num);

    bool single_connection = false;
    if (single_connection)
    {
        if (Kernels.size() != 2)  cout << "more kernel" << endl;
        int neg_num = 0;
        int neg_num2 = 0;
        for (int idx = 0; idx < grid_num; ++idx) {
            Eigen::Vector3d p = GV.row(idx);
            //SDF_gaussian_tubes(idx) = min(SDF_gaussian(idx), generate_tube(p, Kernels[0], Kernels[1], 0.5, 0.5));
            //SDF_gaussian_tubes(idx) = smooth_UnionSDF(SDF_gaussian(idx), generate_tube(p, Kernels[0], Kernels[1], gauss_iso, 0.2), smooth_t);
            //SDF_gaussian_tubes(idx) = generate_tube(p, Kernels[0], Kernels[1], gauss_combined, 0.5);


            //SDF_gaussian_tubes(idx) = generate_tube2(p, Kernels[0], Kernels[1]);
            SDF_gaussian_tubes2(idx) = smooth_UnionSDF(SDF_gaussian(idx), generate_tube2(p, Kernels[0], Kernels[1], gauss_iso, tube_radius), smooth_t);

            /*if (SDF_gaussian_tubes(idx) < isolevel)
                cout << neg_num++ << " sdf1  value: " << SDF_gaussian_tubes(idx) <<endl;*/

        }
    }
    else
    {
        std::cout << "Kernel size: " << Kernels.size() << endl;
		Tube_edges = pores_connection_mst(Kernels, degree_limit);
        for (int idx = 0; idx < grid_num; ++idx) {
            Eigen::Vector3d p = GV.row(idx);
            //SDF_gaussian_tubes(idx) = -max(-SDF_gaussian(idx), generate_tube(p, Kernels[0], Kernels[1], 0.5, 0.5));
            double sdf_p = 1000.0;
            for (auto& e : Tube_edges)
            {
                //sdf_p = min(sdf_p, generate_tube2(p, Kernels[e.from], Kernels[e.to], gauss_iso) );
                sdf_p = min(sdf_p, generate_tube(p, Kernels[e.from], Kernels[e.to], gauss_iso, tube_radius));
            }
            //for (int kx = 0; kx < Kernels.size() - 1; kx++)
            //{
            //    sdf_p = min(sdf_p, generate_tube(p, Kernels[kx], Kernels[kx + 1], gauss_iso, 0.2));
            //    //sdf_p = min(sdf_p, generate_tube2(p, Kernels[kx], Kernels[kx + 1], gauss_iso));
            //}
            SDF_gaussian_tubes2(idx) = smooth_UnionSDF(SDF_gaussian(idx), sdf_p, smooth_t);

            //SDF_gaussian_tubes(idx) = generate_tube(p, Kernels[0], Kernels[1], gauss_iso, 0.5);
        }
    }


	/*string filename1 = "D://VSprojects//TaihuStone//result//sdf1.txt";
	string filename2 = "D://VSprojects//TaihuStone//result//sdf2.txt";
    exportSDF(SDF_gaussian_tubes, filename1);
    exportSDF(SDF_gaussian_tubes2, filename2);*/


    Eigen::MatrixXd V_t; //输出网格顶点
    Eigen::MatrixXi F_t; // 输出网格面片
    MarchingCubes(SDF_gaussian_tubes2, GV, resolution, resolution, resolution, isolevel, V_t, F_t);
    view_model(V_t, F_t);


	// ------------check single component--------------
    //int comp = single_component(V_g, F_g);
    

    //view_model(V_g, F_g);
	//std::string filename = "output/gaussian_pores.stl";
 //   std::filesystem::path filePath(filename);
 //   std::filesystem::path dir = filePath.parent_path();
 //   if (!dir.empty() && !std::filesystem::exists(dir)) {
 //       std::filesystem::create_directories(dir);
 //   }

    // 保存STL文件
    //bool stl_success = igl::writeSTL(filename, V, F, N, igl::FileEncoding::Binary);

    //view_three_models(V_out, F_out, V_g, F_g, V_t, F_t, Eigen::RowVector3d(1, 0, 0));
    view_three_models(V_out, F_out, V_t, F_t, V_t, F_t, Eigen::RowVector3d(1, 0, 0));
    //view_two_models(V_out, F_out, V_g, F_g, Eigen::RowVector3d(1, 0,0));
    // calculate normal
    //Eigen::MatrixXd N;
    //igl::per_face_normals(V_out, F_out, N);
    std::cout << "Generated mesh: " << V_out.rows() << " vertices, " << F_out.rows() << " faces" << std::endl;
    return;


}

void ModelGenerator::sample_interior_points(std::vector<Eigen::Vector3d>& pore_centers, std::vector<double>& pore_sdfs, std::vector<int>& inside_indices, int pores, std::mt19937& gen)
{
    Eigen::VectorXd SDF = this->SDF_ini;
    Eigen::MatrixXd GV = this->GV;
    int grid_num = SDF.size();
    // search inside points
    for (int idx = 0; idx < grid_num; ++idx) {
        if (SDF(idx) < isolevel) {
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

        // 随机选择一个内部点，保存其sdf值
        int chosen_idx = inside_indices[index_dist(gen)];
        Eigen::Vector3d candidate_center = GV.row(chosen_idx).transpose();

        // 检查与已有空洞中心的最小距离
        bool valid = true;
        Eigen::Vector3d min_pt = GV.colwise().minCoeff();
        Eigen::Vector3d max_pt = GV.colwise().maxCoeff();
        Eigen::Vector3d box_size = max_pt - min_pt;
        double volume = box_size.x() * box_size.y() * box_size.z();
        safe_distance = Safe_distance_ratio * std::cbrt(volume / pores);

        for (const auto& existing_center : pore_centers) {
            if ((candidate_center - existing_center).squaredNorm() < safe_distance * safe_distance) {
                valid = false;
                break;
            }
        }

        if (valid) {
            pore_centers.push_back(candidate_center);
            pore_sdfs.push_back(SDF_ini(chosen_idx));
            if(debug_show)
                cout << "pore_sdfs: " << SDF_ini(chosen_idx) << endl;
        }
    }

    std::cout << "Generate " << pore_centers.size() << " kernels" << std::endl;
}

void ModelGenerator::generate_gaussians(std::vector<Eigen::Vector3d> pore_centers, std::vector<double> pore_sdfs, std::mt19937& gen)
{
    Kernels.clear();
    surface_kernels.clear();
    std::uniform_real_distribution<double> dist_amp(amplitude_min, amplitude_max);
    std::uniform_real_distribution<double> dist_sigma(sigma_min, sigma_max);

    // 为每个空洞中心生成随机参数
    int pore_size = pore_centers.size();
    for (size_t i = 0; i < pore_size; ++i) {
        double sigma_val = dist_sigma(gen);
        double amplitude_val = dist_amp(gen);
        //std::cout << "i: " << i << "  " << sigma_val<<"  "<<amplitude_val << std::endl;
        GaussianKernel kernel(pore_centers[i], sigma_val, amplitude_val, pore_sdfs[i]);
        Kernels.push_back(kernel);
        if (kernel.on_surface) surface_kernels.push_back(i);
        /*pore_amplitudes.push_back(dist_amp(gen));
        pore_sigmas.push_back(dist_sigma(gen));*/
    }

    std::cout << "--------------------1. Combine " << Kernels.size() << " Gaussian fileds with " << surface_kernels.size() << " kernels on the surface--------------------" << std::endl;
    if(debug_show)
    {
        cout << "surface index:  ";
        for (auto i : surface_kernels) cout << i << "  ";
        cout << endl;
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
    //std::cout << "total_void: " << total_void << std::endl;
    return  C - total_void;  // 当前使用：负的空洞总和 
}

void ModelGenerator::show_model()
{
    view_two_models(V_ini, F_ini, V_out, F_out, Eigen::RowVector3d(1, 0.0, 0.0));

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
    if(debug_show)
        for (Edge e : mst_edges)
            std::cout << "Edge from " << e.from << " to " << e.to << " with weight " << e.weight << " with length " << e.length << std::endl;
    return mst_edges;
}

std::vector<std::vector<int>> ModelGenerator::construct_adj_list(std::vector<Edge> edges_list, int kernel_num)
{
    std::vector<std::vector<int>> adj_list(kernel_num);
    for (const auto& edge : edges_list) {
        // 由于最小生成树是无向图，一条边代表双向连接。
        // 我们需要将 `to` 添加到 `from` 的邻居列表，同时将 `from` 添加到 `to` 的邻居列表。
        if (edge.from < kernel_num && edge.to < kernel_num) {
            adj_list[edge.from].push_back(edge.to);
            adj_list[edge.to].push_back(edge.from);
        }
    }
    //for (auto adj1 : adj_list)
    //{
    //    cout << "  adj1: " << adj1.size() << " ";
    //    for (auto adj2 : adj1)
    //        cout << "  adj2: " << adj2 << " ";
    //}  

	return adj_list;
}

std::vector<double> ModelGenerator::construct_dist_map(int p_index, std::vector<std::vector<int>> adj)
{
    int n = Kernels.size();
    std::priority_queue<NodeDist, std::vector<NodeDist>, std::greater<NodeDist>> pq;
    std::vector<double> dist_map(n, INF);
    dist_map[p_index] = 0.0;
    pq.push({ p_index, 0.0 });

    while (!pq.empty()) {
        NodeDist current = pq.top();
        pq.pop();

        int u = current.id;
        double d = current.dist;

        if (d > dist_map[u]) continue;

        for (const auto& neighbor : adj[u]) {
            int v = neighbor;
			double length = length_path(u, v);

            if (dist_map[u] + length < dist_map[v]) {
                dist_map[v] = dist_map[u] + length;
                pq.push({ v, dist_map[v] });
            }
        }
    }
	return dist_map;
}


std::vector<int> ModelGenerator::find_path_in_tree(int p1, int p2,  int num_nodes)
{
    // --- 步骤 1: 将边列表转换为邻接表 ，其中 adj_list[i] 将存储所有与节点 i 直接相连的节点的ID。
    if (num_nodes <= 0) {
        return {}; // 如果没有节点，则不可能有路径
    }
    if (p1 == p2) {
        // 处理起点和终点是同一个节点的特殊情况
        return {p1}; // 如果节点重复，则输出p1
    }
    std::vector<std::vector<int>> adj_list = Adj_list;

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
        // 这里体现了邻接表的效率优势：直接访问 adj_list[u] 即可。
        for (int v : adj_list[u]) {
            if (!visited[v]) {
                visited[v] = true; // 标记为已访问
                parent[v] = u;    // 记录父节点，用于后续路径重构
                q.push(v);        // 将邻居节点加入待访问队列
            }
        }
    }

    // --- 步骤 3: 路径重构 ---
    // 如果找到了路径（即BFS到达了终点），我们从终点开始，
    // 利用 `parent` 数组反向回溯，直到回到起点。

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

double ModelGenerator::length_graph_path(int p1, int p2)
{
    std::vector<Edge> mst = Tube_edges;
    if (p1 == p2) return 0.0;
    if (mst.empty()) return -1.0; // 错误：没有树

    std::vector<int> path_ = find_path_in_tree(p1, p2, Kernels.size());

    double length_graph = 0.0;
    for (int i = 0; i < path_.size() - 1; i++) 
    {
        {
            length_graph += distance(Kernels[path_[i]].center, Kernels[path_[i + 1]].center);
        }
        //length_graph += mst[p].length;
        //length_graph += Tube_edges[p].weight;
    }
    return length_graph; // 如果图不连通（理论上MST应该是连通的），返回-1
}

double ModelGenerator::length_path(int p1, int p2)
{
    return distance(Kernels[p1].center, Kernels[p2].center);
}

int  ModelGenerator::find_edge_by_nodes(int from_node, int to_node, const std::vector<Edge> edge_list)
{
    for (int i = 0; i < edge_list.size(); i++)
    {
        if ((edge_list[i].from == from_node && edge_list[i].to == to_node) ||
            (edge_list[i].from == to_node && edge_list[i].to == from_node))
        {
            return i;
        }
    }
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

    double r0 = sqrt(mid_radius_factor) * k1.sigma;
    double r4 = sqrt(mid_radius_factor) * k2.sigma  ;

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


double ModelGenerator::generate_tube2( Eigen::Vector3d& p,   GaussianKernel& k1,   GaussianKernel& k2, double iso_level_C, double mid_radius_factor)
{

    const Eigen::Vector3d& c1 = k1.center;
    const Eigen::Vector3d& c2 = k2.center;

    // 1. 计算 ω_i 和 ω_j, 这里根据sigma进行转换
// 假设 ω = 1/(2 * sigma^2)，保持与标准高斯函数一致
    double omega1 = 1.0 / (2.0 * k1.sigma * k1.sigma);
    double omega2 = 1.0 / (2.0 * k2.sigma * k2.sigma);
    double avg_omega = (omega1 + omega2) / 2.0;
    double mu = 10* mid_radius_factor * avg_omega;  //mid_radius_factor越小，圆越小

    // 2. 计算点到线段的距离（论文中的 ||p - s_ij||）
    Eigen::Vector3d line_dir = c2 - c1;
    double line_length = line_dir.squaredNorm();

    Eigen::Vector3d w = p - c1;
     double dis_c = w.dot(line_dir);
     double dis;
     if (dis_c <= 0) {
         dis = w.squaredNorm();;
     }
     else if (dis_c >= line_dir.dot(line_dir)) {
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
     double tubeValue = tunnelMain + pore1 + pore2;
     //cout << "tunnelMain:  " << tunnelMain << "   " << pore1 << "   " << pore2 << endl;
     //return  tubeValue;
     return iso_level_C - tubeValue;
}

double ModelGenerator::calculate_edge_weight(GaussianKernel k1, GaussianKernel k2)
{
	vector<double> weights{ 0.7, 1.0, 1.3 };  // 分类系数权重：边界-内部，内部-内部，边界-边界
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
    double translucency_score = 1.0;
    int psize = path.size();
    if (psize < 2) {
        cout << "Warnning: illegal path!" << endl;
        return 0.0; // 单个点
    }
    if (psize == 2) 
    {
        //是否横贯模型，如果是返回1，否则返回距离
        int start_ = path[0];
        int end_ = path[1];
        int sam_num = 3;
        double thres = min(Kernels[start_].center_value, Kernels[end_].center_value);
        translucency_score = line_cross_surface(Kernels[start_].center, Kernels[end_].center, thres, sam_num);
    }
    else
    {
        std::vector<GaussianKernel> all_nodes = Kernels;
        int count_inner = 0;
        vector<double> angle_degrees;
        for (size_t i = 1; i < psize - 1; ++i)
        {
            Vector3d prev = all_nodes[path[i - 1]].center;
            Vector3d curr = all_nodes[path[i]].center;
            Vector3d next = all_nodes[path[i + 1]].center;
            if (!all_nodes[path[i]].on_surface) count_inner++;
            double angle_deg = abs_angle(prev - curr, next - curr) / 180.0;
            translucency_score *= angle_deg;
            angle_degrees.push_back(angle_deg);
        }
        if (show_debug)
        {
            cout << "angle_deg: ";
            for (auto ang_ : angle_degrees)
                cout << ang_ << "  ";
            cout << endl;
        }
        translucency_score = std::pow(translucency_score, 1.0 / (psize - 2 + count_inner));// / Kernels.size();
        //translucency_score = std::pow(translucency_score, 1.0 / (psize - 2)) * (1+ 1.0* count_inner / psize);// / Kernels.size();
    }
    
    return translucency_score * log(psize);
}


double ModelGenerator::cal_kernel_translucency(int p_index, int & max_s1, int & max_s2 )  //计算点p所在的所有路径中通透性最大的一条路径作为通透性值
{
    double ave_perm = 0.0;
    int count_ = 0;
    double max_perm = 0.0;
    max_s1 = -1;
    max_s2 = -1;
    std::vector<double> dist_map = construct_dist_map(p_index, Adj_list);
    // 双重循环遍历所有 s1, s2 组合
    // 复杂度 O(K^2)，其中 K 是 surface_points 的数量
    //for (int s1 : surface_kernels)
    for (int i = 0; i < surface_kernels.size(); i++) 
    {
        // 剪枝：如果 s1 无法到达 p_index，则跳过
        int s1 = surface_kernels[i];
        if (dist_map[s1] == INF) continue;
        for (int j = i+1; j < surface_kernels.size(); j++) 
        //for (int s2 : surface_kernels) 
        {
            int s2 = surface_kernels[j];
            // 剪枝：如果 s2 无法到达 p_index，或者 s1==s2 (距离为0)，跳过
            if (dist_map[s2] == INF || s1 == s2) continue;

            // 1. 计算图上距离 (Graph Distance) 路径: s1 -> p_index -> s2
            double graph_dist = dist_map[s1] + dist_map[s2];

            if (graph_dist < 1e-9) continue;

            // 2. 计算欧氏距离 (Euclidean Distance)
            double euclidean_dist = distance(Kernels[s1].center, Kernels[s2].center);

            std::vector<int> path_ = find_specified_path(p_index, s1, s2);
			double path_translucency = calculate_path_translucency(path_);
            if(debug_show)
                cout << "s1: " << s1 << "  s2: " << s2 << endl
                << "angle trans: " << path_translucency << "   length ratio: " << euclidean_dist / graph_dist << endl;
            // 3. 计算通透性
            double perm = path_translucency * euclidean_dist / graph_dist;
            ave_perm += perm;
            count_++;
            // 4. 更新最大值
            if (perm > max_perm)
            {
                max_perm = perm;
                max_s1 = s1;
                max_s2 = s2;
            }
        }
    }

    //show max case
    //std::vector<int> path_ = find_specified_path(p_index, max_s1, max_s2);
    //double path_score = calculate_path_translucency(path_);

    //return ave_perm / count_;
    return max_perm;
}

double ModelGenerator::cal_total_translucency(std::vector<GaussianKernel> gau, std::vector<int> surface_ks)
{
    cout << "--------------------3. Calculating total translucency of mst --------------------" << endl;
    int kernels_num = Kernels.size();
    double total_score = 0.0;
    vector<pair<int, int>> max_paths_kernel;
    //for (int i = 0; i < 1; i++)
    for (int i = 0; i < kernels_num; i++)
    {
        int start = -1, end = -1;
        double score_p = cal_kernel_translucency(i, start, end);
        total_score += score_p;
        max_paths_kernel.push_back(make_pair(start, end));

        cout << "Kernel " << i << " :  max translucency: " << score_p << "   from " << start << " to " << end << endl;
        debug_show = false;
        std::vector<int> path_ = find_specified_path(i, start, end, debug_show);
        double path_translucency = calculate_path_translucency(path_, debug_show);
    }
    return total_score/ kernels_num;
}

std::vector<int> ModelGenerator::find_specified_path(int p_index, int s1, int s2, bool show_debug)
{
    if(s1 ==s2|| s1<0||s2<0) 
    {
        cout << "Warnning: illegal path!" << endl;
        return {};
	}
	int kernel_num = Kernels.size();
    std::vector<int> path1 = find_path_in_tree(p_index, s1, kernel_num);
    //for (auto pp : path1) cout << "path1 : " << pp << "  ";
    std::vector<int> path2 = find_path_in_tree(p_index, s2, kernel_num);
    
    // 合并路径，避免重复包含 p_index
    std::vector<int> full_path(path1.rbegin(), path1.rend());
    full_path.insert(full_path.end(), path2.begin() + 1, path2.end()); // 从 path2 的倒数第二个元素开始添加
    if(show_debug)
    {
        cout << "full path steps: " << full_path.size() << endl;
        for (auto p : full_path) cout << p << " ";
        cout << endl;
    }
    return full_path;
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

    int index = ix * res * res + iy * res + iz;

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
    if (center_sdf + 1e-2 < thres)
    {
        //cout << "the line cross the model!" << endl;
        return 1.0; // 横贯模型
    }
    else
    {
        //cout << "the line on the surafece!" << endl;
        //cout << center_sdf - thres << "  " << abs((center_sdf - thres) / (2 * thres)) << endl;
        return 1.0 - abs((center_sdf - thres) / thres); //表面路径，无内部节点，通透性为1
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
