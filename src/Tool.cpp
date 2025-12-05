#include "Tool.h"


void Mesh2SDF(Eigen::MatrixXd& V,  Eigen::MatrixXi& F, Eigen::MatrixXd& GV, Eigen::VectorXd& SDF)
{
    if (V.rows() == 0 || F.rows() == 0 ) {
        std::cerr << " Input failed! Wrong mesh!" << std::endl;
        return;
    }

    Eigen::VectorXi I;  // 最近面索引
    Eigen::MatrixXd C;  // 最近点坐标
    Eigen::MatrixXd N;  // 内外法向符号

    Eigen::Vector3d bb_min = V.colwise().minCoeff();
    Eigen::Vector3d bb_max = V.colwise().maxCoeff();
    Eigen::Vector3d bb_size = bb_max - bb_min;
    
    //normalize: scaling and moving
    double max_dim = bb_size.maxCoeff();
    double scale_factor = 1.0 / (max_dim + 1e-6);
    Eigen::Vector3d bb_center = (bb_min + bb_max) / 2.0;
    V = (V.rowwise() - bb_center.transpose()) * scale_factor;
    //bb_min = (bb_min - bb_center) * scale_factor;
    bb_min = Vector3d(-0.5, -0.5, -0.5);
    //cout << "bb_min: " << bb_min << endl;
    //bb_min = V.colwise().minCoeff();
    
    int res = Resolution;  // 网格分辨率，可调高以提升精度
    //double dx = bb_size.maxCoeff() / (res - 1);
    double dx = 1.0 / (res - 1);

	int total_points = res * res * res;
    GV.resize(total_points, 3);
    SDF.resize(total_points);
    double void_count = 0;
 
   std:: vector<Eigen::Vector3d> grid_points;
    grid_points.reserve(res * res * res);
    //
    for (int z = 0; z < res; ++z) {
        for (int y = 0; y < res; ++y) {
            for (int x = 0; x < res; ++x) {
                Eigen::Vector3d p = bb_min + Eigen::Vector3d(x, y, z) * dx;
                grid_points.push_back(p);
            }
        }
    }
    //cout << "min: " << grid_points[0] << "    max:  " << grid_points[grid_points.size() - 1] << endl;
    for (int i = 0; i < grid_points.size(); ++i)
        GV.row(i) = grid_points[i];

     std::cout << "Sample resolution: " << res<<"  "<< res<<"  "<< res<<" = "<< GV.rows() << std::endl;

    // 调用 libigl 的 signed_distance()

    igl::signed_distance( GV, V, F, igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER, SDF, I, C, N );

    //Eigen::MatrixXd V_t; //输出网格顶点
    //Eigen::MatrixXi F_t; // 输出网格面片
    //MarchingCubes(SDF, GV, res, res, res, 0, V_t, F_t);  //gaussian combined with tubes
    //view_model(V_t, F_t);

}

bool saveMesh(std::string filename, Eigen::MatrixXd V, Eigen::MatrixXi F)
{
    // 确保输出目录存在
    std::filesystem::path filePath(filename);
    std::filesystem::path dir = filePath.parent_path();
    if (!dir.empty() && !std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
    }
    igl::write_triangle_mesh(filename, V, F);
}

//平滑并集  softmin , k越大，平滑效果越小，趋近于普通并集
double smooth_UnionSDF(double sdf1, double sdf2, double k)   //k y
{
    return -1.0 / k * log(exp(-sdf1 * k) + exp(-sdf2 * k));
}
// SDF平滑交集  softmax, k越大，平滑效果越小，趋近于普通交集
double smooth_IntersecSDF(double sdf1, double sdf2, double k)
{
    return  1.0 / k * log(exp(sdf1 * k) + exp(sdf2 * k));
}

// SDF布尔运算：并集
double unionSDF(double sdf1, double sdf2) 
{
    return std::min(sdf1, sdf2);
}

// SDF布尔运算：交集
double intersectionSDF(double sdf1, double sdf2) 
{
    return std::max(sdf1, sdf2);
}

// SDF布尔运算：差集 (A - B)
double differenceSDF(double sdf1, double sdf2) 
{
    return std::max(sdf1, -sdf2);
}

void MarchingCubes(Eigen::VectorXd& S, Eigen::MatrixXd& GV, int nx, int ny, int nz, double isovalue,  Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    igl::marching_cubes(S, GV, nx, ny, nz, isovalue, V, F);
    if (V.rows() == 0 || F.rows() == 0) {
        std::cerr << "Marching Cubes failed: Empty mesh" << std::endl;
        return;
    }
}


void view_model(Eigen::MatrixXd V1, Eigen::MatrixXi F1)
{
    std::cout << "show libigl viewer" << std::endl;
    igl::opengl::glfw::Viewer viewer;
    //viewer.core().set_uicontrol(false);

    viewer.callback_key_pressed = [](igl::opengl::glfw::Viewer&, unsigned int, int)->bool
        {
            return true; // 阻止默认 key handler 执行，usage 不会打印
        };

    viewer.data().set_mesh(V1, F1);
    viewer.data().show_lines = true;   // 不显示网格线
    //viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色

    viewer.data().point_size = 10; // 让点更显眼
    viewer.launch();
}

void view_two_models(Eigen::MatrixXd V1, Eigen::MatrixXi F1, Eigen::MatrixXd V2, Eigen::MatrixXi F2, Eigen::RowVector3d shift)
{
    std::cout << "show libigl viewer" << std::endl;
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_key_pressed = [](igl::opengl::glfw::Viewer&, unsigned int, int)->bool
        {
            return true; // 阻止默认 key handler 执行，usage 不会打印
        };
    viewer.data().set_mesh(V1, F1);
    viewer.data().show_lines = true;   // 不显示网格线
    //viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色

    int id2 = viewer.append_mesh();
    Eigen::MatrixXd V_shifted = V2;
    V_shifted.rowwise() += shift;  // 向右移动 1 个单位

    viewer.data(id2).set_mesh(V_shifted, F2);
    viewer.data(id2).set_colors(Eigen::RowVector3d(0.8, 0.1, 0.1));

    // 添加辅助点 (高斯核的中心)，设置为红色
   // viewer.data().add_points(kernel_points, Eigen::RowVector3d(1, 0, 0));
    viewer.data().point_size = 10; // 让点更显眼
    viewer.launch();

}

void view_three_models(Eigen::MatrixXd V1, Eigen::MatrixXi F1, Eigen::MatrixXd V2, Eigen::MatrixXi F2, Eigen::MatrixXd V3, Eigen::MatrixXi F3, Eigen::RowVector3d shift)
{
    std::cout << "show libigl viewer" << std::endl;
    igl::opengl::glfw::Viewer viewer;
    viewer.callback_key_pressed = [](igl::opengl::glfw::Viewer&, unsigned int, int)->bool
        {
            return true; // 阻止默认 key handler 执行，usage 不会打印
        };
    viewer.data().set_mesh(V1, F1);
    viewer.data().show_lines = true;   // 不显示网格线
    //viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色

    int id2 = viewer.append_mesh();
    Eigen::MatrixXd V_shifted = V2;
    //V_shifted.rowwise() -= shift;  // 向左移动 1 个单位

    viewer.data(id2).set_mesh(V_shifted, F2);
    viewer.data(id2).set_colors(Eigen::RowVector3d(0.8, 0.1, 0.1));

	int id3 = viewer.append_mesh();
    Eigen::MatrixXd V_shifted3 = V3;
    V_shifted3.rowwise() += shift;  // 向右移动 1 个单位

    viewer.data(id3).set_mesh(V_shifted3, F3);
    viewer.data(id3).set_colors(Eigen::RowVector3d(0.8, 0.8, 0.8));
    // 添加辅助点 (高斯核的中心)，设置为红色
   // viewer.data().add_points(kernel_points, Eigen::RowVector3d(1, 0, 0));
    viewer.data().point_size = 10; // 让点更显眼
    viewer.launch();
}

int single_component(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
    if (F.rows() == 0) return -1;    
    // 构建顶点到面的邻接
    std::vector<std::vector<int>> vertex_faces(V.rows());
    vertex_faces.reserve(V.rows());
    for (int fi = 0; fi < F.rows(); ++fi) {
        for (int k = 0; k < 3; ++k) {
            int v = F(fi, k);
            if (v >= 0 && v < (int)vertex_faces.size())
                vertex_faces[v].push_back(fi);
        }
    }
    // 面的访问标记
    std::vector<char> visited(F.rows(), 0);
    std::vector<int> comp_ids(F.rows(), -1);
    int comp = 0;
    std::vector<int> stack;
    stack.reserve(F.rows());
    std::vector<int> comp_sizes;
    comp_sizes.reserve(16);
    for (int fi = 0; fi < F.rows(); ++fi) {
        if (visited[fi]) continue;
        int size = 0;
        stack.clear();
        stack.push_back(fi);
        visited[fi] = 1;
        comp_ids[fi] = comp;
        while (!stack.empty()) {
            int cur = stack.back(); stack.pop_back();
            ++size;
            // 通过共享顶点扩展
            for (int k = 0; k < 3; ++k) {
                int v = F(cur, k);
                for (int nf : vertex_faces[v]) {
                    if (!visited[nf]) {
                        visited[nf] = 1;
                        comp_ids[nf] = comp;
                        stack.push_back(nf);
                    }
                }
            }
        }
        comp_sizes.push_back(size);
        ++comp;
    }
    if (comp <= 1)
    {
        std::cout << "Only single component" << std::endl;
        return 1; // 已经是单一组件
    }
    else
    {
        std::cout << "Still have " << comp << " components" << std::endl;
        return comp;
    }
     
    // 找最大组件
    //int max_comp = 0;
    //int max_size = comp_sizes[0];
    //for (int i = 1; i < comp; ++i) {
    //    if (comp_sizes[i] > max_size)
    //    {
    //        max_size = comp_sizes[i];
    //        max_comp = i;
    //    }
    //}
    //// 收集保留面
    //std::vector<int> newF_indices;
    //newF_indices.reserve(max_size);
    //for (int fi = 0; fi < F.rows(); ++fi)
    //    if (comp_ids[fi] == max_comp)
    //        newF_indices.push_back(fi);
    //// 复制面
    //Eigen::MatrixXi F_new(newF_indices.size(), 3);
    //for (int i = 0; i < (int)newF_indices.size(); ++i)
    //    F_new.row(i) = F.row(newF_indices[i]);
    //// 压缩顶点：找出现的顶点
    //std::vector<char> used(V.rows(), 0);
    //for (int i = 0; i < F_new.rows(); ++i)
    //    for (int k = 0; k < 3; ++k)
    //        used[F_new(i, k)] = 1;
    //std::vector<int> old2new(V.rows(), -1);
    //int nv = 0;
    //for (int i = 0; i < (int)used.size(); ++i)
    //    if (used[i])
    //        old2new[i] = nv++;
    //Eigen::MatrixXd V_new(nv, 3);
    //for (int i = 0; i < V.rows(); ++i)
    //    if (used[i])
    //        V_new.row(old2new[i]) = V.row(i);
    //for (int i = 0; i < F_new.rows(); ++i)
    //    for (int k = 0; k < 3; ++k)
    //        F_new(i, k) = old2new[F_new(i, k)];
    //V = std::move(V_new);
    //F = std::move(F_new);
    //std::cout << "keepLargestConnectedComponent: kept component " << max_comp << " with " << max_size << " faces." << std::endl;
}

double abs_angle(Vector3d v1, Vector3d v2)
{
    double mag1 = v1.norm();
    double mag2 = v2.norm();
    double angle_deg = 0.0;
    if (mag1 < 1e-9 || mag2 < 1e-9) {
		cout << "Warining: zero length vector in cos_angle!" << endl;
    }
    else
    {
        double dot_product = (v1).dot(v2);
        double cos_theta = dot_product / (mag1 * mag2);
        //cout << "dot_product: " << dot_product << "   cos_theta: " << cos_theta << endl;
        // 夹角可能因浮点误差略超出[-1, 1]范围
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

        double angle_rad = std::acos(cos_theta);
        angle_deg = angle_rad * 180.0 / M_PI;
        //cout << "angle_rad: " << angle_rad << "   angle_deg: " << angle_deg << endl;
    }	
    return angle_deg;
}

double distance(Vector3d v1, Vector3d v2)
{
	return (v1 - v2).norm();
}

double squared_distance(Vector3d v1, Vector3d v2)
{
    return (v1 - v2).squaredNorm();
}

bool align_models_with_pca(const std::string& model1_path, const std::string& model2_path, const std::string& output_path) 
{
    // 读取模型
    MatrixXd V1, V2;
    MatrixXi F1, F2;
    MatrixXd N1, N2;

    if (!igl::read_triangle_mesh(model1_path, V1, F1) ||
        !igl::read_triangle_mesh(model2_path, V2, F2)) {
        return false;
    }
    cout << "Load meshes!" << endl;
    // 计算主成分分析(PCA)获取主方向
    auto compute_principal_axes = [](const MatrixXd& V) -> Matrix3d {
        // 中心化
        Vector3d center = V.colwise().mean();
        MatrixXd centered = V.rowwise() - center.transpose();
        cout << "centering!" << endl;
        // 计算协方差矩阵
        Matrix3d cov = centered.transpose() * centered / (V.rows() - 1);

        // 特征分解
        SelfAdjointEigenSolver<Matrix3d> eigen_solver(cov);
        Matrix3d eigenvectors = eigen_solver.eigenvectors();
        cout << "feature decomposed!" << endl;
        // 按特征值降序排列
        Vector3d eigenvalues = eigen_solver.eigenvalues();
        Matrix3d sorted_eigenvectors;
        for (int i = 0; i < 3; ++i) {
            int max_idx;
            eigenvalues.maxCoeff(&max_idx);
            sorted_eigenvectors.col(2 - i) = eigenvectors.col(max_idx);
            eigenvalues(max_idx) = -1; // 标记为已处理
        }

        return sorted_eigenvectors;
        };

    // 获取主方向
    Matrix3d axes1 = compute_principal_axes(V1);
    Matrix3d axes2 = compute_principal_axes(V2);

    bool flip_x = false, flip_y = false, flip_z = true;

    if (flip_x) axes2.col(0) = -axes2.col(0);
    if (flip_y) axes2.col(1) = -axes2.col(1);
    if (flip_z) axes2.col(2) = -axes2.col(2);

    // 构建旋转矩阵
    Matrix3d rotation = axes1 * axes2.transpose();

    // 先应用旋转
    MatrixXd V2_rotated = V2 * rotation.transpose();
    cout << "Rotated!" << endl;
    // 计算旋转后的包围盒
    Vector3d min1 = V1.colwise().minCoeff();
    Vector3d max1 = V1.colwise().maxCoeff();
    Vector3d min2_rotated = V2_rotated.colwise().minCoeff();
    Vector3d max2_rotated = V2_rotated.colwise().maxCoeff();

    Vector3d center1 = (min1 + max1) / 2.0;
    Vector3d center2_rotated = (min2_rotated + max2_rotated) / 2.0;
    Vector3d size1 = max1 - min1;
    Vector3d size2_rotated = max2_rotated - min2_rotated;

    // 计算缩放和平移
    Vector3d scale = size1.array() / size2_rotated.array();
    Vector3d translation = center1 - center2_rotated;
    cout << "Compute scaling ratio!" << endl;
    // 构建完整变换矩阵
    Matrix4d transform = Matrix4d::Identity();
    transform.block<3, 3>(0, 0) = rotation * scale.asDiagonal();
    transform.block<3, 1>(0, 3) = translation;

    // 应用变换
    MatrixXd V2_transformed(V2.rows(), 3);
    for (int i = 0; i < V2.rows(); ++i) {
        Vector4d v;
        v << V2(i, 0), V2(i, 1), V2(i, 2), 1.0;
        Vector4d v_transformed = transform * v;
        V2_transformed.row(i) = v_transformed.head<3>();
    }

    // 保存结果
    return igl::write_triangle_mesh(output_path, V2_transformed, F2);
}


//B_i,n(t) = C(n, i) * t^i * (1 - t)^(n-i)
double bernstein_basis(int i, int n, double t)
{
    // 组合数 nCk
    auto comb = [](int n, int k)->double {
        if (k < 0 || k > n) return 0.0;
        k = std::min(k, n - k);
        double r = 1.0;
        for (int j = 1; j <= k; ++j) { r *= double(n - (k - j)) / double(j); }
        return r;
        };
    double b = comb(n, i) * std::pow(t, i) * std::pow(1.0 - t, n - i);
    return b;
}

bool exportSDF(Eigen::VectorXd& sdf, std::string& filename) 
{
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Warning: cannot open file! " << filename << std::endl;
        return false;
    }
    outFile << std::fixed << std::setprecision(8);

    for (int i = 0; i < sdf.size(); ++i) 
    {
        outFile << sdf(i);

        // 检查是否是当前行的第100个元素
        if ((i + 1) % 100 == 0) {
            // 如果是，并且不是整个向量的最后一个元素，则换行
            if ((i + 1) != sdf.size()) {
                outFile << '\n';
            }
        }
        else {
            // 如果不是行尾，并且不是整个向量的最后一个元素，则加空格
            if ((i + 1) != sdf.size()) {
                outFile << ' ';
            }
        }
    }
    outFile << '\n';

    return true;
}

void show_path(std::vector<int> path)
{
    cout << "path with "<<path.size()<<" steps : ";
    for (auto p : path)
        cout << p << "  ";
    cout << endl;
}



void geometry_analyzer(Eigen::VectorXd SDF, int resolution, double thres_degree, int overhang_count, int floating_count, std::vector<uint8_t> &overhang_mask, std::vector<uint8_t> &floating_mask)
{
    double overhang_threshold = -std::cos(thres_degree * M_PI / 180.0f);
    overhang_mask.clear(); // 1表示该位置存在悬垂违规
    floating_mask.clear(); // 1表示该位置是悬空孤岛
    overhang_count = 0;
    floating_count = 0;

    int total_voxels = resolution * resolution * resolution;
    overhang_mask.resize(total_voxels, 0);
    floating_mask.resize(total_voxels, 0);

    for (int z = 1; z < resolution - 1; ++z) {
        for (int y = 1; y < resolution - 1; ++y) {
            for (int x = 1; x < resolution - 1; ++x) {
                int idx = x + y * resolution + z* resolution* resolution;
                double val = SDF[idx];

                // 这里的 0.5 是假设体素大小为1，只检测表面附近的体素
                // 在实际应用中，通常检测 abs(val) < voxel_size * sqrt(3)
                if (std::abs(val) < 0.8f) {
                    Vector3d normal = computeGradient(x, y, z, resolution, SDF);

                    // 检查法线Z分量
                    // normal.z < 0 表示朝下
                    // normal.z < -0.707 表示角度陡于45度
                    if (normal.z() < overhang_threshold) {
                        overhang_mask[idx] = 1;
                        overhang_count++;
                    }
                }
            }
        }
    }
    // 2. 检测悬空特征 (Floating Islands)
        // 使用 BFS 从 Z=0 开始标记
    std::vector<bool> visited(total_voxels, false);
    std::queue<int> q;

    // A. 初始化种子：将 Z=0 平面上的所有实体体素加入队列
    for (int y = 0; y < resolution; ++y) {
        for (int x = 0; x < resolution; ++x) {
            int idx = x+ y* resolution;
            // 假设 SDF <= 0 表示实体
            if (SDF[idx] <= 0.0f) {
                q.push(idx);
                visited[idx] = true;
            }
        }
    }

    // B. 广度优先搜索 (BFS) 传播支撑
    int dx[] = { 1, -1, 0, 0, 0, 0 };
    int dy[] = { 0, 0, 1, -1, 0, 0 };
    int dz[] = { 0, 0, 0, 0, 1, -1 };

    while (!q.empty()) {
        int curr_idx = q.front();
        q.pop();

        int cx, cy, cz;
        getCoord(curr_idx, cx, cy, cz);

        for (int i = 0; i < 6; ++i) {
            int nx = cx + dx[i];
            int ny = cy + dy[i];
            int nz = cz + dz[i];

            if (nx>=0&&nx< resolution&& ny >= 0 && ny < resolution&& nz >= 0 && nz < resolution)
            {
                int n_idx = nx+ ny*resolution+ nz* resolution* resolution;
                // 如果邻居是实体 且 未被访问
                if (SDF[n_idx] <= 0.0f && !visited[n_idx]) {
                    visited[n_idx] = true;
                    q.push(n_idx);
                }
            }
        }
    }

    // C. 寻找未被访问的实体 (即悬空部分)
    for (int i = 0; i < total_voxels; ++i) {
        if (SDF[i] <= 0.0f && !visited[i]) {
            floating_mask[i] = 1;
            floating_count++;
        }
    }

}

Vector3d computeGradient(int x, int y, int z, int res, Eigen::VectorXd SDF)
{
    // 边界检查略，循环中已控制范围
    double dx = SDF[x + 1 + y * res + z * res * res] - SDF[x - 1 + y * res + z * res * res];
    double dy = SDF[x + (y + 1) * res + z * res * res] - SDF[x + (y - 1) * res + z * res * res];
    double dz = SDF[x + y * res + (z + 1) * res * res] - SDF[x + y * res + (z - 1) * res * res];

    Vector3d normal(dx, dy, dz);
    double normal_l = normal.norm();
    if(normal_l < 1e-9)
		return Vector3d(0, 0, 0);
    else
        return normal / normal_l;
}

void getCoord(int idx, int res, int& x, int& y, int& z) 
{
    x = idx % res;
    int temp = idx / res;
    y = temp % res;
    z = temp / res;
}

//TO
double smoothHeaviside(double s, double eps)
{
    s = -s;
    if (s > eps) return 1.0;
    if (s < -eps) return 0.0;

    return 0.5 + s / (2.0 * eps)
        + 0.5 / M_PI * sin(M_PI * s / eps);
}

double hardTrans(double s, double iso)
{
    if (s < iso) return 1.0;
    else
        return 0.0;
}

VoxelGrid SDFtoVoxel(std::function<double(const Eigen::Vector3d&)> sdf, Eigen::Vector3d minBox, Eigen::Vector3d maxBox, int nx, int ny, int nz, double eps)
{
    VoxelGrid grid;
    grid.nx = nx; grid.ny = ny; grid.nz = nz;
    grid.origin = minBox;

    grid.dx = (maxBox.x() - minBox.x()) / (nx - 1);
    grid.dy = (maxBox.y() - minBox.y()) / (ny - 1);
    grid.dz = (maxBox.z() - minBox.z()) / (nz - 1);

    grid.rho.resize(nx * ny * nz);

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
                Eigen::Vector3d x(
                    minBox.x() + i * grid.dx,
                    minBox.y() + j * grid.dy,
                    minBox.z() + k * grid.dz
                );

                double phi = sdf(x);   // 你的 SDF 查询
                double rho = smoothHeaviside(phi, eps);
                //double rho = hardTrans(phi, 0.0);

                grid.at(i, j, k) = rho;
            }

    return grid;
}

void saveRawDensity(const VoxelGrid& grid, const std::string& filename)
{
    std::ofstream out(filename, std::ios::binary);
    out.write(reinterpret_cast<const char*>(grid.rho.data()),
        grid.rho.size() * sizeof(double));
    out.close();
}



