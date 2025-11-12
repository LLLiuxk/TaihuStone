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
    for (int i = 0; i < res; ++i) {
        for (int j = 0; j < res; ++j) {
            for (int k = 0; k < res; ++k) {
                Eigen::Vector3d p = bb_min + Eigen::Vector3d(i, j, k) * dx;
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
    F.col(0).swap(F.col(1));  //翻转marching cubes生成的法线方向
}


void view_model(Eigen::MatrixXd V1, Eigen::MatrixXi F1)
{
    std::cout << "show libigl viewer" << std::endl;
    igl::opengl::glfw::Viewer viewer;

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

    // 设置输出精度，以避免科学记数法并保证数据一致性
    outFile << std::fixed << std::setprecision(8);

    for (int i = 0; i < sdf.size(); ++i) {
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
    // 文件末尾可以留一个换行符，这通常是好的实践
    outFile << '\n';

    return true;
}

//visualization
bool compareSDFAndVisualize(
    const Eigen::VectorXd& sdf1,
    const Eigen::VectorXd& sdf2,
    int resX, int resY, int resZ,
    double tolerance,
    const std::string& outputHtmlFile)
{
    // 2. 验证数据
    if (sdf1.size() != sdf2.size()) {
        std::cerr << "Error：two SDF size is not consistent! file1: " << sdf1.size() << ", file2: " << sdf2.size() << std::endl;
        return false;
    }
    long long expected_size = (long long)resX * resY * resZ;
    if (sdf1.size() != expected_size) {
        std::cerr << "错误：SDF数据量 (" << sdf1.size()
            << ") 与提供的分辨率 (" << resX << "x" << resY << "x" << resZ
            << " = " << expected_size << ") 不匹配。" << std::endl;
        return false;
    }

    // 3. 比较并收集误差点
    std::vector<ErrorPoint> errors;
    for (int i = 0; i < sdf1.size(); ++i) {
        double diff = std::abs(sdf1(i) - sdf2(i));
        if (diff > tolerance) {
            ErrorPoint ep;
            ep.value1 = sdf1(i);
            ep.value2 = sdf2(i);
            ep.error = diff;

            // 从1D索引i计算3D坐标(x, y, z)
            // 假设数据是 Z-major order: for z { for y { for x { ... } } }
            ep.z = i / (resX * resY);
            int slice_idx = i % (resX * resY);
            ep.y = slice_idx / resX;
            ep.x = slice_idx % resX;

            errors.push_back(ep);
        }
    }

    // 4. 生成HTML文件
    std::ofstream htmlFile(outputHtmlFile);
    if (!htmlFile.is_open()) {
        std::cerr << "错误：无法创建HTML输出文件: " << outputHtmlFile << std::endl;
        return false;
    }

    // --- 开始写入HTML内容 ---
    htmlFile << R"(<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SDF 差异可视化</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; margin: 0; background-color: #f0f2f5; color: #333; }
        .container { max-width: 1200px; margin: 20px auto; padding: 20px; background-color: #fff; box-shadow: 0 2px 8px rgba(0,0,0,0.1); border-radius: 8px; }
        h1, h2 { color: #1a237e; border-bottom: 2px solid #3949ab; padding-bottom: 10px; }
        canvas { display: block; width: 100%; height: 60vh; background-color: #111; cursor: grab; border-radius: 4px; }
        canvas:active { cursor: grabbing; }
        #tooltip { position: fixed; display: none; background: rgba(255, 255, 255, 0.9); color: #000; padding: 8px 12px; border-radius: 4px; font-size: 14px; pointer-events: none; box-shadow: 0 1px 4px rgba(0,0,0,0.2); }
        table { width: 100%; border-collapse: collapse; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 10px; text-align: left; }
        th { background-color: #e8eaf6; }
        tr:nth-child(even) { background-color: #f9f9f9; }
    </style>
</head>
<body>
    <div class="container">
        <h1>SDF 差异可视化报告</h1>
        <p><strong>网格分辨率:</strong> )" << resX << "x" << resY << "x" << resZ << R"(</p>
        <p><strong>误差阈值:</strong> )" << tolerance << R"(</p>
        <p><strong>发现差异点总数:</strong> )" << errors.size() << R"(</p>
)";

    if (errors.empty()) {
        htmlFile << "<h2>结论：在给定的误差阈值内，两个SDF文件没有差异。</h2>";
    }
    else {
        htmlFile << R"(
        <h2>3D 交互式可视化</h2>
        <p>拖动鼠标以旋转立方体。将鼠标悬停在红点上以查看详细信息。</p>
        <canvas id="sdfCanvas"></canvas>
        <div id="tooltip"></div>

        <h2>差异点详情列表</h2>
        <table>
            <thead>
                <tr>
                    <th>坐标 (X, Y, Z)</th>
                    <th>文件1的值</th>
                    <th>文件2的值</th>
                    <th>绝对误差</th>
                </tr>
            </thead>
            <tbody>
)";
        // 写入表格数据
        for (const auto& err : errors) {
            htmlFile << "<tr><td>(" << err.x << ", " << err.y << ", " << err.z << ")</td><td>"
                << err.value1 << "</td><td>" << err.value2 << "</td><td>" << err.error << "</td></tr>\n";
        }
        htmlFile << R"(
            </tbody>
        </table>
)";
    }

    htmlFile << R"(
    </div>
    <script>
        const canvas = document.getElementById('sdfCanvas');
        const ctx = canvas.getContext('2d');
        const tooltip = document.getElementById('tooltip');

        const resolution = { x: )" << resX << ", y: " << resY << ", z: " << resZ << R"( };
        const errors = [
)";
    // 将误差点数据序列化为JS对象数组
    for (size_t i = 0; i < errors.size(); ++i) {
        const auto& err = errors[i];
        htmlFile << "            { x: " << err.x << ", y: " << err.y << ", z: " << err.z
            << ", val1: " << err.value1 << ", val2: " << err.value2 << ", error: " << err.error << " }";
        if (i < errors.size() - 1) htmlFile << ",\n";
    }
    htmlFile << R"(
        ];

        if (errors.length > 0) {
            let angleX = 0.4;
            let angleY = -0.4;
            let zoom = 1.5;
            let isDragging = false;
            let lastMouseX = 0, lastMouseY = 0;

            const cubeVertices = [
                {x: 0, y: 0, z: 0}, {x: resolution.x, y: 0, z: 0},
                {x: resolution.x, y: resolution.y, z: 0}, {x: 0, y: resolution.y, z: 0},
                {x: 0, y: 0, z: resolution.z}, {x: resolution.x, y: 0, z: resolution.z},
                {x: resolution.x, y: resolution.y, z: resolution.z}, {x: 0, y: resolution.y, z: resolution.z}
            ];

            const cubeEdges = [
                [0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4],
                [0, 4], [1, 5], [2, 6], [3, 7]
            ];

            function project(point3d, w, h) {
                const center = { x: resolution.x / 2, y: resolution.y / 2, z: resolution.z / 2 };
                let x = point3d.x - center.x;
                let y = point3d.y - center.y;
                let z = point3d.z - center.z;

                // Rotate X
                let tempY = y * Math.cos(angleX) - z * Math.sin(angleX);
                let tempZ = y * Math.sin(angleX) + z * Math.cos(angleX);
                y = tempY; z = tempZ;
                
                // Rotate Y
                let tempX = x * Math.cos(angleY) - z * Math.sin(angleY);
                tempZ = x * Math.sin(angleY) + z * Math.cos(angleY);
                x = tempX; z = tempZ;

                const scale = Math.min(w, h) / Math.max(resolution.x, resolution.y, resolution.z) * 0.8 * zoom;
                
                return {
                    x: w / 2 + x * scale,
                    y: h / 2 + y * scale,
                    z: z, // For depth sorting/sizing
                    scale: scale
                };
            }

            function draw() {
                const w = canvas.width;
                const h = canvas.height;
                ctx.clearRect(0, 0, w, h);

                // Draw cube edges
                ctx.strokeStyle = 'rgba(100, 150, 255, 0.5)';
                ctx.lineWidth = 1;
                cubeEdges.forEach(edge => {
                    const p1 = project(cubeVertices[edge[0]], w, h);
                    const p2 = project(cubeVertices[edge[1]], w, h);
                    ctx.beginPath();
                    ctx.moveTo(p1.x, p1.y);
                    ctx.lineTo(p2.x, p2.y);
                    ctx.stroke();
                });

                // Draw error points
                const projectedErrors = errors.map(p => ({ ...p, proj: project(p, w, h) }));
                projectedErrors.sort((a, b) => b.proj.z - a.proj.z); // Depth sort

                projectedErrors.forEach(p => {
                    ctx.beginPath();
                    ctx.arc(p.proj.x, p.proj.y, 3, 0, 2 * Math.PI);
                    ctx.fillStyle = 'rgba(255, 50, 50, 0.8)';
                    ctx.fill();
                });
            }

            function resizeCanvas() {
                canvas.width = canvas.clientWidth;
                canvas.height = canvas.clientHeight;
                draw();
            }

            canvas.addEventListener('mousedown', e => {
                isDragging = true;
                lastMouseX = e.clientX;
                lastMouseY = e.clientY;
            });
            window.addEventListener('mouseup', () => isDragging = false);
            window.addEventListener('mousemove', e => {
                if (!isDragging) return;
                const dx = e.clientX - lastMouseX;
                const dy = e.clientY - lastMouseY;
                angleY += dx * 0.01;
                angleX -= dy * 0.01;
                lastMouseX = e.clientX;
                lastMouseY = e.clientY;
                draw();
            });

            canvas.addEventListener('wheel', e => {
                e.preventDefault();
                zoom *= (1 - e.deltaY * 0.001);
                zoom = Math.max(0.1, Math.min(10, zoom));
                draw();
            });

            canvas.addEventListener('mousemove', e => {
                if (isDragging) return;
                const rect = canvas.getBoundingClientRect();
                const mouseX = e.clientX - rect.left;
                const mouseY = e.clientY - rect.top;

                let hoveredPoint = null;
                const w = canvas.width;
                const h = canvas.height;
                
                const projectedErrors = errors.map(p => ({ ...p, proj: project(p, w, h) }));
                projectedErrors.sort((a, b) => b.proj.z - a.proj.z);

                for (const p of projectedErrors) {
                    const dist = Math.sqrt((p.proj.x - mouseX)**2 + (p.proj.y - mouseY)**2);
                    if (dist < 8) {
                        hoveredPoint = p;
                        break;
                    }
                }

                if (hoveredPoint) {
                    tooltip.style.display = 'block';
                    tooltip.style.left = `${e.clientX + 15}px`;
                    tooltip.style.top = `${e.clientY}px`;
                    tooltip.innerHTML = `
                        <strong>Coord:</strong> (${hoveredPoint.x}, ${hoveredPoint.y}, ${hoveredPoint.z})<br>
                        <strong>File 1:</strong> ${hoveredPoint.val1.toFixed(6)}<br>
                        <strong>File 2:</strong> ${hoveredPoint.val2.toFixed(6)}<br>
                        <strong>Error:</strong> ${hoveredPoint.error.toFixed(6)}
                    `;
                } else {
                    tooltip.style.display = 'none';
                }
            });

            window.addEventListener('resize', resizeCanvas);
            resizeCanvas();
        }
    </script>
</body>
</html>
)";
    // --- HTML内容写入结束 ---

    htmlFile.close();
    std::cout << "可视化报告已生成: " << outputHtmlFile << std::endl;
    return true;
}

