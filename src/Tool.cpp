#include "Tool.h"

vector<RowVector3d> colors = { RowVector3d(0.90, 0.62, 0.0), //orange
    RowVector3d(0.9, 0.75, 0.40), //light orange
    RowVector3d(0.2, 0.44, 0.90),  //blue
    RowVector3d(0.14, 0.24, 0.52),  //drak blue
    RowVector3d(0.8, 0.33, 0.20),  //muted red
    RowVector3d(0.8, 0.1, 0.1),  //red
    RowVector3d(0.75, 0.75, 0.75),  // dark gray
    RowVector3d(0.3, 0.68, 0.45),  //drak green
    RowVector3d(0.0, 0.62, 0.45),  //drak green
    RowVector3d(0.8, 0.47, 0.65), //pink red
    RowVector3d(0.84, 0.37, 0.0), //brick red
    RowVector3d(0.94, 0.89, 0.26),  //yellow
    RowVector3d(0.90, 0.80, 0.60)  //yellow
};



double Mesh2SDF(Eigen::MatrixXd& V,  Eigen::MatrixXi& F, Eigen::MatrixXd& GV, Eigen::VectorXd& SDF, Eigen::Vector3d& bb_min, Eigen::Vector3d& bb_max)
{
    if (V.rows() == 0 || F.rows() == 0 ) {
        std::cerr << " Input failed! Wrong mesh!" << std::endl;
        return 0.0;
    }

    Eigen::VectorXi I;  // 最近面索引
    Eigen::MatrixXd C;  // 最近点坐标
    Eigen::MatrixXd N;  // 内外法向符号

    bb_min = V.colwise().minCoeff();
    bb_max = V.colwise().maxCoeff();
    //cout << bb_min << endl << bb_max << endl;
    Eigen::Vector3d bb_size = bb_max - bb_min;
    
    //normalize: scaling and moving
    double max_dim = bb_size.maxCoeff();
    double scale_factor = 1.0 / (max_dim + 1e-6);
    Eigen::Vector3d bb_center = (bb_min + bb_max) / 2.0;
    V = (V.rowwise() - bb_center.transpose()) * scale_factor;
    bb_min = (bb_min - bb_center) * scale_factor;
    bb_max = (bb_max - bb_center) * scale_factor;
    Vector3d left_corner(-0.5, -0.5, -0.5);

    std::string ori_model = "D:/VSprojects/TaihuStone/model/ori/" + input_file + ".stl";
    saveMesh(ori_model, V, F);
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
                Eigen::Vector3d p = left_corner + Eigen::Vector3d(x, y, z) * dx;
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
    return scale_factor;
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
    cout << "Save meth to " << filename << "  successfully" << endl;
	return true;
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


void view_model(Eigen::MatrixXd V1, Eigen::MatrixXi F1, std::string win_name, bool show_line)
{
    std::cout << "show libigl viewer" << std::endl;
    igl::opengl::glfw::Viewer viewer;

    viewer.core().background_color.setConstant(1.0f); // White background

    viewer.data().set_mesh(V1, F1);
    viewer.data().show_lines = show_line;   // 不显示网格线
    viewer.data().show_faces = !show_line;   // 不显示三角面
    //viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色

    viewer.data().point_size = 10; // 让点更显眼
    //viewer.launch();

    //viewer.core().align_camera_center(V1, F1);
    float deg2rad = float(M_PI) / 180.0f;
    //show_degree_x = -90, show_degree_z =128
    Eigen::Quaternionf q =
        Eigen::AngleAxisf(show_degree_x * deg2rad, Eigen::Vector3f::UnitX()) *
        Eigen::AngleAxisf(show_degree_y * deg2rad, Eigen::Vector3f::UnitY()) *
        Eigen::AngleAxisf(show_degree_z * deg2rad, Eigen::Vector3f::UnitZ());
    viewer.core().trackball_angle = q;


    viewer.launch(false, win_name);
}

void view_two_models(Eigen::MatrixXd V1, Eigen::MatrixXi F1, Eigen::MatrixXd V2, Eigen::MatrixXi F2, Eigen::RowVector3d shift, std::string win_name)
{
    std::cout << "show libigl viewer" << std::endl;
    igl::opengl::glfw::Viewer viewer;
    viewer.core().background_color.setConstant(1.0f); // White background

    viewer.data().set_mesh(V1, F1);
    //viewer.data().show_lines = false;   // 不显示网格线
    viewer.data().show_lines = true;   // 不显示网格线
    viewer.data().show_faces = false;   // 不显示三角面
    viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.8, 0.8));
    viewer.data().shininess = 500.0;
    //viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色

    int id2 = viewer.append_mesh();
    viewer.data(id2).set_mesh(V2, F2);
    //viewer.data(id3).set_colors(Eigen::RowVector3d(0.4, 0.4, 0.2));
    viewer.data(id2).show_lines = false;   // 不显示网格线

    //viewer.data(id2).set_colors(Eigen::RowVector3d(0.95, 0.7, 0.30));
    viewer.data(id2).set_colors(Eigen::RowVector3d(0.18, 0.67, 0.97));
    viewer.data(id2).shininess = 2000.0;
    // 添加辅助点 (高斯核的中心)，设置为红色
   // viewer.data().add_points(kernel_points, Eigen::RowVector3d(1, 0, 0));
    //viewer.data().point_size = 10; // 让点更显眼

    float deg2rad = float(M_PI) / 180.0f;
    Eigen::Quaternionf q =
        Eigen::AngleAxisf(show_degree_x * deg2rad, Eigen::Vector3f::UnitX()) *
        Eigen::AngleAxisf(show_degree_y * deg2rad, Eigen::Vector3f::UnitY()) *
        Eigen::AngleAxisf(show_degree_z * deg2rad, Eigen::Vector3f::UnitZ());
    viewer.core().trackball_angle = q;

    viewer.launch(false, win_name);
}

void view_three_models(Eigen::MatrixXd V1, Eigen::MatrixXi F1, Eigen::MatrixXd V2, Eigen::MatrixXi F2, Eigen::MatrixXd V3, Eigen::MatrixXi F3, Eigen::RowVector3d shift, std::string win_name)
{
    std::cout << "show libigl viewer" << std::endl;
    igl::opengl::glfw::Viewer viewer;
    viewer.core().background_color.setConstant(1.0f); // White background

    viewer.data().set_mesh(V1, F1);
    viewer.data().show_lines = false;   // 不显示网格线
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
    viewer.launch(false, win_name);
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

Eigen::Matrix3d construct_R(double u, double v, double theta_max_rad, Eigen::Vector3d z_axis)
{
    double cos_theta = 1.0 - u * (1.0 - std::cos(theta_max_rad));
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));

    Eigen::Vector3d d(
        sin_theta * std::cos(v),
        sin_theta * std::sin(v),
        cos_theta
    );
    d.normalize();
    // 如果几乎没有旋转，直接返回单位矩阵
    if ((d - z_axis).norm() < 1e-6) {
        return Eigen::Matrix3d::Identity();
    }
    // Rodrigues 公式
    Eigen::Vector3d axis = z_axis.cross(d);
    axis.normalize();

    double angle = std::acos(
        std::min(1.0, std::max(-1.0, z_axis.dot(d)))
    );

    Eigen::Matrix3d K;
    K << 0, -axis.z(), axis.y(),
        axis.z(), 0, -axis.x(),
        -axis.y(), axis.x(), 0;

    Eigen::Matrix3d R =
        Eigen::Matrix3d::Identity()
        + std::sin(angle) * K
        + (1.0 - std::cos(angle)) * K * K;

    return R;
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

void vis_Kernels_Tubes(std::vector<Eigen::Vector3d>& points, std::vector<pair<int, int>>& connections, std::string win_name)
{
    const double sphereRadius = 0.01;
    const double lineWidth = 7;
    const int  sphereSubdiv = 3;

    // ---- Color palette ----
    const RowVector3d sphereColor = colors[12]; // muted red
    const RowVector3d lineColor= colors[6]; // dark gray
    //const RowVector3d sphereColor(0.75, 0.33, 0.20); // muted red
    //const RowVector3d lineColor(0.20, 0.20, 0.20); // dark gray

    igl::opengl::glfw::Viewer viewer;
	viewer.core().background_color.setConstant(1.0f); // White background

    MatrixXd V_unit;
    MatrixXi F_unit;
    igl::icosahedron(V_unit, F_unit);
    for (int i = 0; i < sphereSubdiv; ++i)
        igl::loop(V_unit, F_unit, V_unit, F_unit);

    for (int i = 0; i < V_unit.rows(); ++i)
        V_unit.row(i).normalize(); // project to unit sphere
    V_unit *= sphereRadius;

    int Nv = V_unit.rows();
    int Nf = F_unit.rows();
    int Np = points.size();
    MatrixXd V_all(Np * Nv, 3);
    MatrixXi F_all(Np * Nf, 3);

    for (int i = 0; i < Np; ++i)
    {
        V_all.block(i * Nv, 0, Nv, 3) =
            V_unit.rowwise() + points[i].transpose();

        F_all.block(i * Nf, 0, 
            Nf, 3) =
            F_unit.array() + i * Nv;
    }
 //   std::string outputPrefix = "D:/VSprojects/TaihuStone/result/" + input_file + "_" + std::to_string(PoresNum) + "_opt/";
	//saveMesh(outputPrefix + "kernels.stl", V_all, F_all);
    viewer.data().set_mesh(V_all, F_all);
    viewer.data().set_colors(sphereColor);
    viewer.data().show_lines = false;
    viewer.data().shininess = 200.0;

    // Draw lines based on the connections
    Eigen::MatrixXd P1(connections.size(), 3);
    Eigen::MatrixXd P2(connections.size(), 3);
    for (int i = 0; i < connections.size(); ++i) {
        P1.row(i) = points[connections[i].first].transpose();
        P2.row(i) = points[connections[i].second].transpose();
    }
    viewer.data().add_edges(P1, P2, lineColor);

    // Set line width
    viewer.data().line_width = lineWidth;

    //viewer.core().align_camera_center(V1, F1);

    // ===精确设置视角===
    float deg2rad = float(M_PI) / 180.0f;
    Eigen::Quaternionf q =
        Eigen::AngleAxisf(show_degree_x * deg2rad, Eigen::Vector3f::UnitX()) *
        Eigen::AngleAxisf(show_degree_y * deg2rad, Eigen::Vector3f::UnitY()) *
        Eigen::AngleAxisf(show_degree_z * deg2rad, Eigen::Vector3f::UnitZ());
    viewer.core().trackball_angle = q;

    // Launch the viewer
    viewer.launch(false, win_name + to_string(connections.size()));
}

void vis_compare_cons(std::vector<Eigen::Vector3d>& points, std::vector<pair<int, int>>& connections, std::vector<int> mask, std::string win_name)
{
    double sphereRadius = 0.01;
    double lineWidth = 7;
    int  sphereSubdiv = 3;

    // ---- Color palette ----
    const RowVector3d sphereColor = colors[12]; // muted red
    const RowVector3d lineColor = colors[6]; // 
    const RowVector3d lineColor1 = colors[5]; // 
    const RowVector3d lineColor2 = colors[7]; // 

    //const RowVector3d sphereColor(0.75, 0.33, 0.20); // muted red
    //const RowVector3d lineColor(0.20, 0.20, 0.20); // dark gray

    igl::opengl::glfw::Viewer viewer;
    viewer.core().background_color.setConstant(1.0f); // White background

    MatrixXd V_unit;
    MatrixXi F_unit;
    igl::icosahedron(V_unit, F_unit);
    for (int i = 0; i < sphereSubdiv; ++i)
        igl::loop(V_unit, F_unit, V_unit, F_unit);

    for (int i = 0; i < V_unit.rows(); ++i)
        V_unit.row(i).normalize(); // project to unit sphere
    V_unit *= sphereRadius;

    int Nv = V_unit.rows();
    int Nf = F_unit.rows();
    int Np = points.size();
    MatrixXd V_all(Np * Nv, 3);
    MatrixXi F_all(Np * Nf, 3);

    for (int i = 0; i < Np; ++i)
    {
        V_all.block(i * Nv, 0, Nv, 3) =
            V_unit.rowwise() + points[i].transpose();

        F_all.block(i * Nf, 0,
            Nf, 3) =
            F_unit.array() + i * Nv;
    }

    viewer.data().set_mesh(V_all, F_all);
    viewer.data().set_colors(sphereColor);
    viewer.data().show_lines = false;
    viewer.data().shininess = 200.0;

    // Draw lines based on the connections
    Eigen::MatrixXd P1(connections.size(), 3);
    Eigen::MatrixXd P2(connections.size(), 3);
    Eigen::MatrixXd C(connections.size(), 3);
    for (int i = 0; i < connections.size(); ++i) {
        P1.row(i) = points[connections[i].first].transpose();
        P2.row(i) = points[connections[i].second].transpose();

        if (mask[i] == 0) {
            C.row(i) = lineColor;
        }
        else if (mask[i] == 1) {
            C.row(i) = lineColor1;
        }
        else if (mask[i] == 2)
        {
            C.row(i) = lineColor2;
        }
     }
    viewer.data().add_edges(P1, P2, C);

    // Set line width
    viewer.data().line_width = lineWidth;

    //viewer.core().align_camera_center(V1, F1);

    // ===精确设置视角===
    float deg2rad = float(M_PI) / 180.0f;
    Eigen::Quaternionf q =
        Eigen::AngleAxisf(show_degree_x * deg2rad, Eigen::Vector3f::UnitX()) *
        Eigen::AngleAxisf(show_degree_y * deg2rad, Eigen::Vector3f::UnitY()) *
        Eigen::AngleAxisf(show_degree_z * deg2rad, Eigen::Vector3f::UnitZ());
    viewer.core().trackball_angle = q;

    // Launch the viewer
    viewer.launch(false, win_name + to_string(connections.size()));
}

static RowVector3d hsv_to_rgb(double h, double s, double v)
{
    h = std::fmod(std::max(0.0, h), 1.0);
    s = std::max(0.0, std::min(1.0, s));
    v = std::max(0.0, std::min(1.0, v));

    double r = v, g = v, b = v;
    if (s > 1e-9) {
        double h6 = h * 6.0;
        int sector = static_cast<int>(std::floor(h6)) % 6;
        double f = h6 - std::floor(h6);
        double p = v * (1.0 - s);
        double q = v * (1.0 - s * f);
        double t = v * (1.0 - s * (1.0 - f));

        switch (sector) {
        case 0: r = v; g = t; b = p; break;
        case 1: r = q; g = v; b = p; break;
        case 2: r = p; g = v; b = t; break;
        case 3: r = p; g = q; b = v; break;
        case 4: r = t; g = p; b = v; break;
        default: r = v; g = p; b = q; break;
        }
    }
    return RowVector3d(r, g, b);
}

static RowVector3d kernel_path_color(int kernel_idx, int kernel_num)
{
    if (kernel_num <= 0) {
        return RowVector3d(0.95, 0.89, 0.26);
    }
    // Golden-ratio hue stepping gives stable, well-spread colors for many kernels.
    constexpr double golden_ratio = 0.6180339887498949;
    double h = std::fmod((kernel_idx * golden_ratio), 1.0);
    return hsv_to_rgb(h, 0.75, 0.95);
}

void vis_translucency_cons(
    std::vector<Eigen::Vector3d>& points,
    std::vector<pair<int, int>>& connections,
    std::vector<int> edge_usage,
    std::string win_name,
    const std::vector<std::vector<int>>& kernel_paths)
{
    double sphereRadius = 0.01;
    double lineWidth = 7;
    int  sphereSubdiv = 3;

    const RowVector3d sphereColor = colors[12];
    const RowVector3d lineColor = colors[6];
    const RowVector3d lineColorUsed = colors[12]; // yellow

    if (edge_usage.size() != connections.size()) {
        cout << "Warning: edge_usage size(" << edge_usage.size()
            << ") != connections size(" << connections.size() << "), fallback to gray lines." << endl;
        edge_usage.assign(connections.size(), 0);
    }

    igl::opengl::glfw::Viewer viewer;
    viewer.core().background_color.setConstant(1.0f);

    MatrixXd V_unit;
    MatrixXi F_unit;
    igl::icosahedron(V_unit, F_unit);
    for (int i = 0; i < sphereSubdiv; ++i)
        igl::loop(V_unit, F_unit, V_unit, F_unit);

    for (int i = 0; i < V_unit.rows(); ++i)
        V_unit.row(i).normalize();
    V_unit *= sphereRadius;

    int Nv = V_unit.rows();
    int Nf = F_unit.rows();
    int Np = points.size();
    MatrixXd V_all(Np * Nv, 3);
    MatrixXi F_all(Np * Nf, 3);

    for (int i = 0; i < Np; ++i)
    {
        V_all.block(i * Nv, 0, Nv, 3) =
            V_unit.rowwise() + points[i].transpose();

        F_all.block(i * Nf, 0, Nf, 3) =
            F_unit.array() + i * Nv;
    }

    viewer.data().set_mesh(V_all, F_all);
    viewer.data().set_colors(sphereColor);
    viewer.data().show_lines = false;
    viewer.data().shininess = 200.0;

    if (!kernel_paths.empty()) {
        // Draw only edges unused by any max-translucency path in gray to avoid
        // complete overlap with colored paths.
        int unused_edge_num = 0;
        for (int i = 0; i < static_cast<int>(connections.size()); ++i) {
            if (i < static_cast<int>(edge_usage.size()) && edge_usage[i] <= 0) {
                unused_edge_num++;
            }
        }

        if (unused_edge_num > 0) {
            Eigen::MatrixXd P1_base(unused_edge_num, 3);
            Eigen::MatrixXd P2_base(unused_edge_num, 3);
            int write_idx = 0;
            for (int i = 0; i < static_cast<int>(connections.size()); ++i) {
                if (i >= static_cast<int>(edge_usage.size()) || edge_usage[i] > 0) {
                    continue;
                }
                P1_base.row(write_idx) = points[connections[i].first].transpose();
                P2_base.row(write_idx) = points[connections[i].second].transpose();
                write_idx++;
            }
            viewer.data().add_edges(P1_base, P2_base, lineColor);
        }

        int valid_kernel_path_num = 0;
        for (size_t k = 0; k < kernel_paths.size(); ++k) {
            const std::vector<int>& path = kernel_paths[k];
            if (path.size() < 2) continue;

            std::vector<Eigen::Vector3d> p1_vec;
            std::vector<Eigen::Vector3d> p2_vec;
            p1_vec.reserve(path.size() - 1);
            p2_vec.reserve(path.size() - 1);

            for (size_t i = 0; i + 1 < path.size(); ++i) {
                int u = path[i];
                int v = path[i + 1];
                if (u < 0 || v < 0 || u >= points.size() || v >= points.size()) {
                    continue;
                }
                p1_vec.push_back(points[u]);
                p2_vec.push_back(points[v]);
            }

            if (p1_vec.empty()) continue;

            ++valid_kernel_path_num;
            Eigen::MatrixXd P1_path(p1_vec.size(), 3);
            Eigen::MatrixXd P2_path(p2_vec.size(), 3);
            for (int i = 0; i < static_cast<int>(p1_vec.size()); ++i) {
                P1_path.row(i) = p1_vec[i].transpose();
                P2_path.row(i) = p2_vec[i].transpose();
            }

            RowVector3d color_k = kernel_path_color(static_cast<int>(k), static_cast<int>(kernel_paths.size()));
            int path_layer = viewer.append_mesh();
            viewer.data(path_layer).add_edges(P1_path, P2_path, color_k);
            viewer.data(path_layer).line_width = lineWidth + 1.5;
        }

        cout << "[vis_translucency_cons] draw " << valid_kernel_path_num
            << " kernel max-translucency paths with unique colors." << endl;
    }
    else {
        // Fallback to original binary coloring mode.
        Eigen::MatrixXd P1(connections.size(), 3);
        Eigen::MatrixXd P2(connections.size(), 3);
        Eigen::MatrixXd C(connections.size(), 3);
        for (int i = 0; i < connections.size(); ++i) {
            P1.row(i) = points[connections[i].first].transpose();
            P2.row(i) = points[connections[i].second].transpose();
            C.row(i) = (edge_usage[i] > 0) ? lineColorUsed : lineColor;
        }
        viewer.data().add_edges(P1, P2, C);
    }
    viewer.data().line_width = lineWidth;

    float deg2rad = float(M_PI) / 180.0f;
    Eigen::Quaternionf q =
        Eigen::AngleAxisf(show_degree_x * deg2rad, Eigen::Vector3f::UnitX()) *
        Eigen::AngleAxisf(show_degree_y * deg2rad, Eigen::Vector3f::UnitY()) *
        Eigen::AngleAxisf(show_degree_z * deg2rad, Eigen::Vector3f::UnitZ());
    viewer.core().trackball_angle = q;

    viewer.launch(false, win_name + to_string(connections.size()));
}

void vis_opt_cons(std::vector<Eigen::Vector3d>& points, std::vector<pair<int, int>>& connections, std::vector<int> mask, std::string win_name)
{
    int delete_in = 0;  //选择显示第几条删除的边
    double sphereRadius = 0.01;
    double lineWidth = 7;
    int  sphereSubdiv = 3;

    // ---- Color palette ----
    const RowVector3d sphereColor = colors[12]; // muted red
    const RowVector3d sphereColor2 = colors[10]; // muted red
    const RowVector3d lineColor = colors[6]; //colors[1]; //  orange: origin
    const RowVector3d lineColor1 = colors[5]; // red: replace -
    const RowVector3d lineColor2 = colors[2]; // blue : replace +
    const RowVector3d lineColor3 = colors[7]; // green : add
    //const RowVector3d sphereColor(0.75, 0.33, 0.20); // muted red
    //const RowVector3d lineColor(0.20, 0.20, 0.20); // dark gray

    igl::opengl::glfw::Viewer viewer;
    viewer.core().background_color.setConstant(1.0f); // White background

    MatrixXd V_unit;
    MatrixXi F_unit;
    igl::icosahedron(V_unit, F_unit);
    for (int i = 0; i < sphereSubdiv; ++i)
        igl::loop(V_unit, F_unit, V_unit, F_unit);

    for (int i = 0; i < V_unit.rows(); ++i)
        V_unit.row(i).normalize(); // project to unit sphere
    V_unit *= sphereRadius;

    int Nv = V_unit.rows();
    int Nf = F_unit.rows();
    int Np = points.size();
    MatrixXd V_all(Np * Nv, 3);
    MatrixXi F_all(Np * Nf, 3);

    for (int i = 0; i < Np; ++i)
    {
        V_all.block(i * Nv, 0, Nv, 3) =
            V_unit.rowwise() + points[i].transpose();

        F_all.block(i * Nf, 0,
            Nf, 3) =
            F_unit.array() + i * Nv;
    }

    // Draw lines based on the connections
    int show_num = PoresNum - 1; // connections.size()
    Eigen::MatrixXd P1(show_num, 3);
    Eigen::MatrixXd P2(show_num, 3);
    //Eigen::MatrixXd C(show_num, 3);
    vector<int> deleted_index; // = 0;
    for (int i = 0; i < show_num; i++) {
        if (mask[i] == 1)
        {
            deleted_index.push_back(i);
        }
        P1.row(i) = points[connections[i].first].transpose();
        P2.row(i) = points[connections[i].second].transpose();
    }

    //找到要绘制的点
    int draw_one = 0;

    if (delete_in >= deleted_index.size()) {
        cout << "beyond delete index!" << endl;
        return;
    }
    //去掉删除的边
    P1.row(deleted_index[delete_in]) = points[connections[deleted_index[delete_in]].first].transpose();
    P2.row(deleted_index[delete_in]) = points[connections[deleted_index[delete_in]].first].transpose();

    pair<int, int> e1 = connections[deleted_index[delete_in]];
    pair<int, int> e2 = connections[show_num + delete_in];

    if(e1.first == e2.first || e1.first == e2.second)
        draw_one = e1.first;
    else
        draw_one = e1.second;

    //cout << "e1: " << deleted_index[delete_in] << "   e2: " << show_num + delete_in << "   draw_one: "<< draw_one<< endl;
    int highlight_id = draw_one;
    MatrixXd V_high = V_unit*1.5;
    V_high.rowwise() += points[highlight_id].transpose();

    viewer.data().set_mesh(V_high, F_unit);
    viewer.data().set_colors(sphereColor2); // 红色
    viewer.data().show_lines = false;
    viewer.data().shininess = 2000.0;
    //viewer.core().align_camera_center(V_all, F_all);
    viewer.data().add_edges(P1, P2, lineColor);

    // Set line width
    viewer.data().line_width = lineWidth;

    vector<int> draw_index;
    draw_index.push_back(deleted_index[delete_in]);
    for (int i = show_num; i < connections.size(); i++)
    {
        //cout << i<<" connections: " << connections[i].first << "   " << connections[i].second << "   mask i :"<<mask[i]<<endl;
        if(connections[i].first == draw_one || connections[i].second == draw_one)
			draw_index.push_back(i);
    }
    int ano_num = draw_index.size();
    Eigen::MatrixXd P3(ano_num, 3);
    Eigen::MatrixXd P4(ano_num, 3);
    Eigen::MatrixXd C(ano_num, 3);
    for (int i = 0; i < ano_num; i++) {
        //if (mask[i] == 3) continue;
        //cout << "draw_index.size: "<<ano_num << "  connections[draw_index[i]] " << connections[draw_index[i]].first << "   " << connections[draw_index[i]].second << "    mask:"<<mask[draw_index[i]]<<endl;
        P3.row(i) = points[connections[draw_index[i]].first].transpose();
        P4.row(i) = points[connections[draw_index[i]].second].transpose();

        if (mask[draw_index[i]] == 1) {
            C.row(i) = lineColor1;
        }
        else if (mask[draw_index[i]] == 2) {
            C.row(i) = lineColor2;
        }
        else 
        {
            C.row(i) = lineColor3;
        }
    }
    
    int id2 = viewer.append_mesh();
    viewer.data(id2).add_edges(P3, P4, C);

    // Set line width
    viewer.data(id2).line_width = 2* lineWidth;

    viewer.data(id2).set_mesh(V_all, F_all);
    viewer.data(id2).set_colors(sphereColor);
    viewer.data(id2).show_lines = false;
    viewer.data(id2).shininess = 200.0;


    //viewer.core().align_camera_center(V1, F1);

    // ===精确设置视角===
    float deg2rad = float(M_PI) / 180.0f;
    Eigen::Quaternionf q =
        Eigen::AngleAxisf(show_degree_x * deg2rad, Eigen::Vector3f::UnitX()) *
        Eigen::AngleAxisf(show_degree_y * deg2rad, Eigen::Vector3f::UnitY()) *
        Eigen::AngleAxisf(show_degree_z * deg2rad, Eigen::Vector3f::UnitZ());
    viewer.core().trackball_angle = q;

    // Launch the viewer
    viewer.launch(false, win_name + to_string(connections.size()));

    //int show_num = 59; // connections.size()
    //Eigen::MatrixXd P1(show_num, 3);
    //Eigen::MatrixXd P2(show_num, 3);
    //Eigen::MatrixXd C(show_num, 3);
    //for (int i = 0; i < show_num; i++) {
    //    //if (mask[i] == 3) continue;
    //    P1.row(i) = points[connections[i].first].transpose();
    //    P2.row(i) = points[connections[i].second].transpose();

    //    if (mask[i] == 0) {
    //        C.row(i) = lineColor;
    //    }
    //    else if (mask[i] == 1) {
    //        C.row(i) = lineColor1;
    //    }
    //    else if (mask[i] == 2)
    //    {
    //        C.row(i) = lineColor2;
    //    }
    //    else if (mask[i] == 3)
    //    {
    //        C.row(i) = lineColor3;
    //    }
    //}
    //viewer.data().add_edges(P1, P2, C);

    //// Set line width
    //viewer.data().line_width = lineWidth;
}


void vis_KerLine_model(Eigen::MatrixXd V1, Eigen::MatrixXi F1, std::vector<Eigen::Vector3d>& points, std::vector<pair<int, int>>& connections, bool show_line, std::string win_name)
{
    const double sphereRadius = 0.01;
    const double lineWidth = 5;
    const int  sphereSubdiv = 3;
    // ---- Color palette ----
    const RowVector3d sphereColor = colors[5]; // muted red
    const RowVector3d lineColor = colors[1]; // dark gray
    //const RowVector3d sphereColor(0.75, 0.33, 0.20); // muted red
    //const RowVector3d lineColor(0.20, 0.20, 0.20); // dark gray

    igl::opengl::glfw::Viewer viewer;
    viewer.core().background_color.setConstant(1.0f); // White background

    MatrixXd V_unit;
    MatrixXi F_unit;
    igl::icosahedron(V_unit, F_unit);
    for (int i = 0; i < sphereSubdiv; ++i)
        igl::loop(V_unit, F_unit, V_unit, F_unit);

    for (int i = 0; i < V_unit.rows(); ++i)
        V_unit.row(i).normalize(); // project to unit sphere
    V_unit *= sphereRadius;

    int Nv = V_unit.rows();
    int Nf = F_unit.rows();
    int Np = points.size();
    MatrixXd V_all(Np * Nv, 3);
    MatrixXi F_all(Np * Nf, 3);

    for (int i = 0; i < Np; ++i)
    {
        V_all.block(i * Nv, 0, Nv, 3) =
            V_unit.rowwise() + points[i].transpose();

        F_all.block(i * Nf, 0,
            Nf, 3) =
            F_unit.array() + i * Nv;
    }

    viewer.data().set_mesh(V_all, F_all);
    viewer.data().set_colors(sphereColor);
    viewer.data().show_lines = false;
    viewer.data().shininess = 200.0;

    // Draw lines based on the connections
    Eigen::MatrixXd P1(connections.size(), 3);
    Eigen::MatrixXd P2(connections.size(), 3);
    for (int i = 0; i < connections.size(); ++i) {
        P1.row(i) = points[connections[i].first].transpose();
        P2.row(i) = points[connections[i].second].transpose();
    }
    viewer.data().add_edges(P1, P2, lineColor);
    viewer.data().line_width = lineWidth;

    // Set line width
    int id2 = viewer.append_mesh();

    viewer.data(id2).set_mesh(V1, F1);
    viewer.data(id2).show_lines = show_line;   // 不显示网格线
    viewer.data(id2).show_faces = !show_line;   // 不显示三角面
	viewer.data(id2).line_color = Eigen::Matrix<float, 4, 1>(0.5f, 0.5f, 0.5f, 1.0f);
    //viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色
    //viewer.data(id2).set_colors(sphereColor);//(Eigen::RowVector3d(0.95, 0.95, 0.95));
    viewer.data(id2).shininess = 500.0;

    //viewer.data(id2).point_size = 10; // 让点更显眼
    //viewer.launch();

    //viewer.core().align_camera_center(V1, F1);
    float deg2rad = float(M_PI) / 180.0f;
    Eigen::Quaternionf q =
        Eigen::AngleAxisf(show_degree_x * deg2rad, Eigen::Vector3f::UnitX()) *
        Eigen::AngleAxisf(show_degree_y * deg2rad, Eigen::Vector3f::UnitY()) *
        Eigen::AngleAxisf(show_degree_z * deg2rad, Eigen::Vector3f::UnitZ());
    viewer.core().trackball_angle = q;

    // Launch the viewer
    viewer.launch(false, win_name + to_string(connections.size()));
}

void geometry_analyzer(Eigen::VectorXd SDF, int resolution, double thres_degree, int &overhang_count, int& floating_count, std::vector<uint8_t> &overhang_mask, std::vector<uint8_t> &floating_mask)
{
    double overhang_threshold = -std::cos(thres_degree * M_PI / 180.0f);
    overhang_mask.clear(); // 1表示该位置存在悬垂违规
    floating_mask.clear(); // 1表示该位置是悬空孤岛
    overhang_count = 0;
    floating_count = 0;

    int total_voxels = resolution * resolution * resolution;
    overhang_mask.resize(total_voxels, 0);
    floating_mask.resize(total_voxels, 0);

    int cal_num = 0;
    for (int z = 1; z < resolution - 1; ++z) {
        for (int y = 1; y < resolution - 1; ++y) {
            for (int x = 1; x < resolution - 1; ++x) {
                int idx = x + y * resolution + z* resolution* resolution;
                double val = SDF[idx];
                cal_num++;
                //if(cal_num%
                // ==0)
                if (val<0)
                    cout << "idx: " << idx << "  val:"<< val<<endl;
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
    cout << "cal_num: " << cal_num << endl;
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
        getCoord(curr_idx, resolution, cx, cy, cz);

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
    cout << "寻找未被访问的实体 (即悬空部分)" << endl;
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


// pca point cloud
Eigen::Vector3d computePrincipalDirection(const std::vector<Eigen::Vector3d>& points)
{
    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    for (auto& p : points) mean += p;
    mean /= points.size();

    Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
    for (auto& p : points) {
        Eigen::Vector3d q = p - mean;
        C += q * q.transpose();
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(C);
    return solver.eigenvectors().col(2).normalized(); // 最大特征值
}

//TO
double smoothHeaviside(double s, double eps)   
{
    if (s < -eps) return 1.0;
    if (s > eps) return 0.0;

    double x = s / eps;
    return 0.5 * (1.0 - sin(M_PI * x / 2.0));
}

double hardTrans(double s, double iso)
{
    if (s < iso) return 1.0;
    else
        return 0.0;
}

VoxelGrid SDFtoVoxel(Eigen::VectorXd& sdf, Eigen::Vector3d minBox, Eigen::Vector3d maxBox, int nx, int ny, int nz)
{
    //cout << "SDFtoVoxel: " << sdf.size()<< endl;
    VoxelGrid grid;
    grid.nx = nx; grid.ny = ny; grid.nz = nz;
    grid.origin = minBox;

    grid.dx = (maxBox.x() - minBox.x()) / (nx - 1);
    grid.dy = (maxBox.y() - minBox.y()) / (ny - 1);
    grid.dz = (maxBox.z() - minBox.z()) / (nz - 1);

    grid.rho.resize(nx * ny * nz);

    int neg_num = 0;
    double eps = 0.2;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
            {
				int index = i + j * nx + k * nx * ny;
                double phi = sdf(index);   // 你的 SDF 查询
                //double rho = smoothHeaviside(phi, eps);
                double rho = hardTrans(phi, 0.0);
                if (rho > 0.9) neg_num++;
                grid.at(i, j, k) = rho;
            }
    //cout << "SDFtoVoxel over!" <<  "   neg_num: " << neg_num << endl; 
    return grid;
}

void saveVoxelToRaw(std::string filename, VoxelGrid& grid)
{
    std::ofstream out(filename, std::ios::binary);
    out.write(reinterpret_cast<const char*>(grid.rho.data()),
        grid.rho.size() * sizeof(double));
    out.close();
}

void saveVoxelToVTK(std::string filename, VoxelGrid& grid)
{
    std::ofstream out(filename);
    out << "# vtk DataFile Version 3.0\n";
    out << "Voxel Density\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS "
        << grid.nx << " "
        << grid.ny << " "
        << grid.nz << "\n";

    out << "ORIGIN "
        << grid.origin.x() << " "
        << grid.origin.y() << " "
        << grid.origin.z() << "\n";

    out << "SPACING "
        << grid.dx << " "
        << grid.dy << " "
        << grid.dz << "\n";

    out << "POINT_DATA " << grid.nx * grid.ny * grid.nz << "\n";
    out << "SCALARS density double\n";
    out << "LOOKUP_TABLE default\n";

    for (double v : grid.rho)
        out << v << "\n";

    out.close();
}

void  saveSDFtoVTI(const std::string filename, Eigen::VectorXd& sdf, int nx, int ny, int nz) {
    std::ofstream out(filename);
    if (!out) return;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <ImageData WholeExtent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    out << "    <Piece Extent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\">\n";
    out << "      <PointData Scalars=\"Density\">\n";
    out << "        <DataArray type=\"Float32\" Name=\"Density\" format=\"ascii\">\n";

	Eigen::VectorXd rho_tilde = sdf;
    
	int n_vars = nx * ny * nz;
    for (int i = 0; i < n_vars; ++i) {
        float rho = smoothHeaviside(rho_tilde(i), 0.1);
        out << rho << " ";
        if (i % 20 == 0) out << "\n";
        /*out << (float)rho_tilde(i) << " ";
        if (i % 20 == 0) out << "\n";*/
    }

    out << "\n        </DataArray>\n";
    out << "      </PointData>\n";
    out << "      <CellData>\n";
    out << "      </CellData>\n";
    out << "    </Piece>\n";
    out << "  </ImageData>\n";
    out << "</VTKFile>\n";

    std::cout << "Saved VTI: " << filename << std::endl;
}

void saveSDFtoNPY(std::string filename, Eigen::VectorXd& sdf, int res) {

    std::vector<double> voxel_grid;
    for (auto sdf_ : sdf)
        voxel_grid.push_back(sdf_);

    saveVoxelGridAsNPY(voxel_grid, res, filename);

}

// 保存体素网格为NPY格式
void saveVoxelGridAsNPY(std::vector<double>& voxel_grid, int res, std::string filename) 
{
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "无法创建NPY文件: " << filename << std::endl;
        return;
    }
    // NPY文件头格式
    // Magic number (6 bytes): \x93NUMPY
    file.write("\x93NUMPY", 6);
    // Version (2 bytes): 1.0
    file.write("\x01\x00", 2);

    // Header dictionary
    std::string dtype = "'<f8'";     // 修改为：小端双精度浮点数
    std::string fortran_order = "False";
    std::string shape = "(" + std::to_string(res) + ", " + std::to_string(res) + ", " + std::to_string(res) + ")";

    std::string header_dict = "{'descr': " + dtype + ", 'fortran_order': " + fortran_order + ", 'shape': " + shape + ", }";

    // 补齐到16字节对齐
    while ((header_dict.length() + 1) % 16 != 15) {
        header_dict += " ";
    }
    header_dict += "\n";


    // Header length (2 bytes, little endian)
    int16_t header_len = static_cast<uint16_t>(header_dict.length());
    file.write(reinterpret_cast<const char*>(&header_len), 2);

    // Header dictionary
    file.write(header_dict.c_str(), header_dict.length());

    std::vector<double> voxel_grid_reordered;
    for (int z = 0; z < res; ++z) {
        for (int y = 0; y < res; ++y) {
            for (int x = 0; x < res; ++x) {
                // 将数据按 z, y, x 顺序排列
                int index = (z * res * res) + (y * res) + x;
                voxel_grid_reordered.push_back(voxel_grid[index]);
            }
        }
    }

    // Data (按照x,y,z顺序存储)
    file.write(reinterpret_cast<const char*>(voxel_grid_reordered.data()), voxel_grid_reordered.size() * sizeof(double)); 

    file.close();
    std::cout << "NPY文件保存成功: " << filename << " (大小: " << res << "x" << res << "x" << res << ")" << std::endl;
}

bool readNPYtoVoxel(const std::string& filename, std::vector<double>& voxel_grid, int& res)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Could not load NPY input: " << filename << std::endl;
        return false;
    }

    // --- 1. 读取并校验 Magic ---
    char magic[6];
    file.read(magic, 6);
    if (std::string(magic, 6) != "\x93NUMPY") return false;

    // --- 2. 跳过 Version (2 bytes) ---
    file.ignore(2);

    // --- 3. 读取 Header 长度 ---
    uint16_t header_len;
    file.read(reinterpret_cast<char*>(&header_len), 2);

    // --- 4. 读取 Header 内容 ---
    std::string header_dict(header_len, ' ');
    file.read(&header_dict[0], header_len);

    // --- 5. 解析 Header (关键修改) ---
    // 检查是 float (<f4) 还是 double (<f8)
    bool is_float32 = false;
    if (header_dict.find("'<f4'") != std::string::npos || header_dict.find("\"<f4\"") != std::string::npos) {
        is_float32 = true;
        //std::cout << "检测到数据格式: Float32 (<f4)" << std::endl;
    }
    else if (header_dict.find("'<f8'") != std::string::npos || header_dict.find("\"<f8\"") != std::string::npos) {
        is_float32 = false;
        //std::cout << "检测到数据格式: Float64/Double (<f8)" << std::endl;
    }
    else {
        //std::cerr << "错误: 未知或不支持的数据类型。" << std::endl;
        return false;
    }

    // 解析 shape
    std::regex shape_regex(R"('shape':\s*\((\d+),\s*(\d+),\s*(\d+)\))");
    std::smatch match;
    if (std::regex_search(header_dict, match, shape_regex)) {
        res = std::stoi(match[1].str()); // 假设是立方体，取第一维
    }
    else {
        std::cerr << "Warnning: could not analyze shape" << std::endl;
        return false;
    }

    size_t total_elements = static_cast<size_t>(res) * res * res;
    voxel_grid.resize(total_elements);

    // --- 6. 分情况读取数据 (关键修改) ---
    if (is_float32) {
        // 如果是 float，先创建一个临时的 float buffer
        std::vector<float> temp_buffer(total_elements);

        file.read(reinterpret_cast<char*>(temp_buffer.data()), total_elements * sizeof(float));

        if (!file) {
            std::cerr << "错误: 文件数据读取不完整 (期望 " << total_elements * sizeof(float) << " 字节)" << std::endl;
            return false;
        }

        // 转换 float -> double
        for (size_t i = 0; i < total_elements; ++i) {
            voxel_grid[i] = static_cast<double>(temp_buffer[i]);
        }
    }
    else {
        // 如果本来就是 double，直接读
        file.read(reinterpret_cast<char*>(voxel_grid.data()), total_elements * sizeof(double));

        if (!file) {
            std::cerr << "错误: 文件数据读取不完整 (期望 " << total_elements * sizeof(double) << " 字节)" << std::endl;
            return false;
        }
    }

    file.close();
    std::cout << "Loading NPY success(Res: " << res << ") from" << filename<<std::endl;
    return true;
}

void writeAdjacencyListToFile(const std::vector<std::vector<int>>& adjList,
    const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    for (const auto& neighbors : adjList) {
        // 写入每个顶点的邻接列表
        for (size_t i = 0; i < neighbors.size(); ++i) {
            outFile << neighbors[i];
            if (i != neighbors.size() - 1) {
                outFile << " ";  // 用空格分隔相邻节点
            }
        }
        outFile << std::endl;  // 每行代表一个顶点的邻接表
    }

    outFile.close();
    std::cout << "邻接表已成功写入到文件: " << filename << std::endl;
}

std::vector<std::vector<int>> readAdjacencyListFromFile(const std::string& filename) {
    std::vector<std::vector<int>> adjList2;
    std::ifstream inFile(filename);

    if (!inFile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return adjList2;
    }

    std::string line;
    while (std::getline(inFile, line)) {
        std::vector<int> neighbors;
        std::stringstream ss(line);
        int value;

        // 从每一行读取所有整数
        while (ss >> value) {
            neighbors.push_back(value);
        }

        adjList2.push_back(neighbors);
    }

    inFile.close();
    return adjList2;
}

//string
std::string to_string_pre(double value, int precision)
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

//double add_iso_surface_noise(
//    const Eigen::Vector3d& p,
//    double sdf_value,
//    double band_width,
//    double noise_amplitude,
//    double spatial_frequency)
//{
//    init_field_noise();
//
//    // 1. 等值面权重
//    double w = std::exp(-(sdf_value * sdf_value) / (2.0 * band_width * band_width)); 
//
//    if (w < 1e-6)
//        return sdf_value;
//
//    // 2. 连续空间噪声
//    double n = g_field_noise.GetNoise(
//        float(p.x() * spatial_frequency),
//        float(p.y() * spatial_frequency),
//        float(p.z() * spatial_frequency)
//    );
//
//    // 3. SDF 扰动
//    return sdf_value + noise_amplitude * w * n;
//}

void add_noise_near_isosurface(
    Eigen::VectorXd& S,              // 标量场（in-place 修改）
    const Eigen::MatrixXd& GV,        // 体素坐标
    double iso_value,                 // MC 等值面
    double noise_amplitude,            // 噪声幅值（建议 0.05 ~ 0.3 * 场尺度）
    double band_width,                 // 等值面带宽
    double spatial_frequency           // 噪声空间频率
)
{
    //init_noise();
    init_field_noise();

    const int N = static_cast<int>(S.size());
    assert(GV.rows() == N);

    for (int i = 0; i < N; ++i)
    {
        double F = S(i);

        // ---- 等值面权重（只在 iso 附近生效）----
        double w = std::exp(
            -std::pow((F - iso_value) / band_width, 2.0)
        );

        if (w < 1e-4)
            continue;

        // ---- 空间噪声 ----
        const Eigen::Vector3d& p = GV.row(i);
        double n = g_field_noise.GetNoise(
            float(p.x() * spatial_frequency),
            float(p.y() * spatial_frequency),
            float(p.z() * spatial_frequency)
        );

        // ---- 标量场扰动（值扰动）----
        S(i) += noise_amplitude * w * n;
    }
}


//cal sigma
double sigma_value(double v, double n, double w, double iso)
{
    double k = 4.0 / 3 * M_PI * pow(-2 * log(iso), 1.5);
    return cbrt(v / (k * n * w));
}

vector<vector<int>> compare_edges(const vector<pair<int, int>>& ini, const vector<pair<int, int>> & final)
{
    vector<vector<int>> mask;
    vector<int> mask1(ini.size(), 1);   // 默认全部标记为删除
    vector<int> mask2(final.size(), 2); // 默认全部标记为新增

    // 统计 final 中各边出现的次数
    map<pair<int, int>, vector<int>> final_map;
    for (int i = 0; i < final.size(); ++i) {
        final_map[final[i]].push_back(i);
    }

    // 遍历 ini，寻找匹配项
    for (int i = 0; i < ini.size(); ++i) {
        auto it = final_map.find(ini[i]);
        if (it != final_map.end() && !it->second.empty()) {
            // 找到匹配：说明这条边“没动”
            mask1[i] = 0;

            // 取出 final 中的索引并标记为 0
            int final_idx = it->second.back();
            mask2[final_idx] = 0;

            // 消耗掉这个索引（处理重复边的情况）
            it->second.pop_back();
        }
    }
    mask.push_back(mask1);
    mask.push_back(mask2);
    int del_num = 0, add_num = 0;
    for (auto m1 : mask1)
        if (m1 == 1) del_num++;
    for (auto m2 : mask2)
        if (m2 == 2) add_num++;
    cout << "delete num : " << del_num << "   add num:" << add_num << endl;

    return mask;
}

vector<int> compare_edges2(const vector<pair<int, int>>& ini, const vector<pair<int, int>> & final, vector<pair<int, int>>& edge_con_total, std::vector<int> rep_vec)
{
    // ini- changed-add: ini顺序与changed顺序不匹配
    edge_con_total.clear();
    vector<int> mask(ini.size(), 0);
    edge_con_total = ini;
    //调整顺序
    std::sort(rep_vec.begin(), rep_vec.end());
    for (auto rv : rep_vec)
    {
		mask[rv] = 1;
		edge_con_total.push_back(final[rv]);
        mask.push_back(2);
    }
    for (int i = ini.size(); i < final.size(); i++)
    {
		edge_con_total.push_back(final[i]);
        mask.push_back(3);
    }

	//cout << "edge_con_total size:" << edge_con_total.size() << "   "<<mask.size()<<endl;
 //   for(int i = 0; i < edge_con_total.size(); i++)
	//	cout << edge_con_total[i].first << "  " << edge_con_total[i].second << "    "<<mask[i]<<endl;
 //   for (int i = 0; i < final.size(); i++)
 //   {
 //       cout << "num " << i << ": " << final[i].first << "  " << final[i].second << "                       ";
 //       if (i < ini.size())
 //           cout << ini[i].first << "  " << ini[i].second << endl;
 //       else
 //           cout << "default" << endl;
 //   }
    return mask;

}



vector<int> cal_max_degree(std::vector<std::vector<int>> Adj_list)
{
    int max_degree = 0;
    int max_degree_num = 0;
    for (const auto& neighbors : Adj_list) {
        int degree = neighbors.size();
        if (degree > max_degree) {
            max_degree = degree;
            max_degree_num = 1;
        }
        else if (degree == max_degree)
            max_degree_num++;
    }
    return { max_degree , max_degree_num };
}

int save_translucency_summary(
    const std::string& output_prefix,
    const std::string& input_name,
    int max_degree,
    int max_degree_num,
    double final_translucency,
    double sum_t_angle,
    double sum_t_length,
    double sum_t_location,
    double sum_s_horiz,
    double sum_weighted,
    int kernel_num,
    int valid_paths,
    const std::vector<double>& kt_weights,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<int>& edge_usage
)
{
    double kernel_num_safe = (kernel_num <= 0) ? 1.0 : static_cast<double>(kernel_num);
    double avg_t_angle = sum_t_angle / kernel_num_safe;
    double avg_t_length = sum_t_length / kernel_num_safe;
    double avg_t_location = sum_t_location / kernel_num_safe;
    double avg_s_horiz = sum_s_horiz / kernel_num_safe;
    double avg_weighted = sum_weighted / kernel_num_safe;

    int used_edge_num = 0;
    for (auto use_times : edge_usage) {
        if (use_times > 0) used_edge_num++;
    }

    cout << "MaxDegree: " << max_degree << "   num: " << max_degree_num << endl;
    cout << "Final translucency score: " << final_translucency << endl;
    cout << "Para sensitivity on max VP paths: "
        << "sum(T_angle)=" << avg_t_angle
        << ", sum(T_length)=" << avg_t_length
        << ", sum(T_location)=" << avg_t_location
        << ", sum(S_horiz)=" << avg_s_horiz
        << ", weighted_sum=" << avg_weighted
        << ", valid_paths=" << valid_paths << endl;
    cout << "KT_weights: ";
    for (size_t i = 0; i < kt_weights.size(); ++i) {
        cout << kt_weights[i];
        if (i + 1 < kt_weights.size()) cout << ", ";
    }
    cout << endl;
    cout << "[Final translucency edges] used " << used_edge_num
        << " / " << edge_usage.size() << " edges." << endl;

    std::filesystem::create_directories(output_prefix);
    std::string trans_txt = output_prefix + input_name + "_translucency_summary.txt";
    std::ofstream trans_ofs(trans_txt, std::ios::out);
    if (!trans_ofs.is_open()) {
        cout << "Warning: failed to open translucency summary file: " << trans_txt << endl;
        return used_edge_num;
    }

    trans_ofs << std::fixed << std::setprecision(6);
    trans_ofs << "MaxDegree: " << max_degree << "   num: " << max_degree_num << "\n";
    trans_ofs << "Final translucency score: " << final_translucency << "\n";
    trans_ofs << "Para sensitivity on max VP paths: "
        << "sum(T_angle)=" << avg_t_angle
        << ", sum(T_length)=" << avg_t_length
        << ", sum(T_location)=" << avg_t_location
        << ", sum(S_horiz)=" << avg_s_horiz
        << ", weighted_sum=" << avg_weighted
        << ", valid_paths=" << valid_paths << "\n";
    trans_ofs << "KT_weights: ";
    for (size_t i = 0; i < kt_weights.size(); ++i) {
        trans_ofs << kt_weights[i];
        if (i + 1 < kt_weights.size()) trans_ofs << ", ";
    }
    trans_ofs << "\n";
    trans_ofs << "[Final translucency edges] used " << used_edge_num
        << " / " << edge_usage.size() << " edges.\n";
    trans_ofs << "Used edges detail (edge_index, from, to, usage_count):\n";

    size_t row_num = std::min(edges.size(), edge_usage.size());
    for (size_t ei = 0; ei < row_num; ++ei) {
        if (edge_usage[ei] > 0) {
            trans_ofs << ei << ", "
                << edges[ei].first << ", "
                << edges[ei].second << ", "
                << edge_usage[ei] << "\n";
        }
    }
    trans_ofs.close();
    return used_edge_num;
}

void append_translucency_summary_metrics(
    const std::string& output_prefix,
    const std::string& input_name,
    int floating_voxel_count,
    int unsupported_voxel_count
)
{
    std::string trans_txt = output_prefix + input_name + "_translucency_summary.txt";
    std::ofstream trans_ofs(trans_txt, std::ios::app);
    if (!trans_ofs.is_open()) {
        cout << "Warning: failed to append translucency summary file: " << trans_txt << endl;
        return;
    }

    trans_ofs << "Floating voxel count (removed from SDF): " << floating_voxel_count << "\n";
    trans_ofs << "Unsupported voxel count: " << unsupported_voxel_count << "\n";
    trans_ofs.close();
}


