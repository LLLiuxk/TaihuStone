//#include <igl/readOFF.h>
//#include <igl/readOBJ.h>
//#include <igl/readSTL.h>
//#include <igl/marching_cubes.h>
//#include <igl/opengl/glfw/Viewer.h>
//#include <igl/copyleft/cgal/mesh_boolean.h>
//#include <igl/read_triangle_mesh.h>
//#include <igl/write_triangle_mesh.h>
//#include <Eigen/Dense>
//#include <iostream>
//#include <vector>
//#include <random>
//
//#include "Tool.h"
//
//// 定义高斯核的数据结构
//struct GaussianKernel {
//    Eigen::Vector3d center; // 核的中心位置
//    double sigma;         // 核的大小/影响力范围 (高斯函数的标准差)
//};
//
//int main(int argc, char* argv[]) {
//    // =========================================================================
//    // 1. 参数设置 (您可以调整这些值来观察不同的效果)
//    // =========================================================================
//    const int num_kernels = 1;      // 您希望在空间中生成的随机核的数量
//    const int grid_resolution = 64; // 采样网格的精细度。越高越精细，但计算越慢。
//    const double isolevel = 0.8;     // 等值面的阈值。决定了"高斯球"的表面在哪里。
//    // 可以尝试 0.5, 1.0 等值来观察变化。
//
//// =========================================================================
//// 2. 生成随机的高斯核 (Kernel_set)
//// =========================================================================
//    std::vector<GaussianKernel> Kernel_set;
//    Kernel_set.reserve(num_kernels);
//
//    // 使用C++11的随机数生成器，比rand()效果更好
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<> pos_dist(-0.8, 0.8); // 位置分布在[-0.8, 0.8]的立方体内
//    std::uniform_real_distribution<> sigma_dist(0.2, 0.4); // 大小(sigma)分布在[0.2, 0.4]之间
//
//    Eigen::MatrixXd kernel_points(num_kernels, 3); // 用于可视化核的中心点
//
//    std::cout << "生成 " << num_kernels << " 个随机高斯核..." << std::endl;
//    for (int i = 0; i < num_kernels; ++i) {
//        GaussianKernel kernel;
//        kernel.center = Eigen::Vector3d(pos_dist(gen), pos_dist(gen), pos_dist(gen));
//        kernel.sigma = sigma_dist(gen);
//        Kernel_set.push_back(kernel);
//        kernel_points.row(i) = kernel.center;
//    }
//
//    // =========================================================================
//    // 3. 创建采样网格和标量场
//    // =========================================================================
//    std::cout << "创建 " << grid_resolution << "x" << grid_resolution << "x" << grid_resolution << " 的采样网格..." << std::endl;
//
//    const int grid_size = pow(grid_resolution, 3);
//    // GV: 存储网格顶点的坐标 (grid_resolution^3 个点)
//    Eigen::MatrixXd GV(grid_size, 3);
//    // S: 存储在每个网格顶点上计算出的标量值
//    Eigen::VectorXd S(grid_size);
//
//
//    // 定义网格的边界，确保能包围所有核
//    double bound = 1.5;
//    Eigen::RowVector3d grid_min(-bound, -bound, -bound);
//    Eigen::RowVector3d grid_max(bound, bound, bound);
//    Eigen::RowVector3d grid_step = (grid_max - grid_min) / (grid_resolution - 1);
//
//    std::cout << "计算标量场 (所有核在每个网格点上的影响力总和)..." << std::endl;
//    for (int k = 0; k < grid_resolution; ++k) {
//        for (int j = 0; j < grid_resolution; ++j) {
//            for (int i = 0; i < grid_resolution; ++i) {
//                // 计算当前网格点的索引
//                int index = i + j * grid_resolution + k * grid_resolution * grid_resolution;
//
//                // 计算当前网格点的世界坐标
//                Eigen::Vector3d p = (grid_min + Eigen::RowVector3d(i * grid_step.x(), j * grid_step.y(), k * grid_step.z())).transpose();
//                GV.row(index) = p;
//
//                // 计算所有高斯核在该点p的影响力总和
//                double total_influence = 0.0;
//                for (const auto& kernel : Kernel_set) {
//                    double sq_dist = (p - kernel.center).squaredNorm();
//                    double sigma_sq = kernel.sigma * kernel.sigma;
//                    // 高斯函数: exp(-d^2 / (2*sigma^2))
//                    total_influence += exp(-sq_dist / (2.0 * sigma_sq));
//                }
//                S(index) = total_influence;
//            }
//        }
//    }
//
//    // =========================================================================
//    // 4. 运行 Marching Cubes 算法生成网格
//    // =========================================================================
//    std::cout << "运行 Marching Cubes 算法 (isolevel = " << isolevel << ")..." << std::endl;
//    Eigen::MatrixXd V; // 输出的网格顶点
//    Eigen::MatrixXi F; // 输出的网格面片
//
//    igl::marching_cubes(S, GV, grid_resolution, grid_resolution, grid_resolution, isolevel, V, F);
//
//    if (V.rows() == 0) {
//        std::cerr << "警告: Marching Cubes未能生成任何顶点。可能是isolevel设置得太高或太低，或者核太分散了。请尝试调整参数。" << std::endl;
//    }
//    else {
//        std::cout << "成功生成网格: " << V.rows() << " 个顶点, " << F.rows() << " 个面片。" << std::endl;
//    }
//
//
//    // --- 1. 加载模型A (您想要从中减去高斯球的模型) ---
//    Eigen::MatrixXd VA;
//    Eigen::MatrixXi FA;
//    // 请将 "path/to/your/model.obj" 替换为您自己的模型文件路径
//    // 这个模型必须是封闭的（水密的）！
//    	
////    if (!igl::readOBJ("D:/VSprojects/TaihuStone/high_rock.obj", VA, FA)) {
//    if (!igl::readOBJ("D:/VSprojects/TaihuStone/high_rock2.obj", VA, FA)) {
//        std::cerr << "Error: Could not load model A." << std::endl;
//        return 1;
//    }
//    std::cout << "Model A loaded successfully." << std::endl;
//
//    Eigen::MatrixXd GV2;
//    Eigen::VectorXd SDF2;
//    OBJ2SDF(VA,FA, GV2, SDF2 );
//    Eigen::MatrixXd V3;
//    Eigen::MatrixXi F3;
//    igl::marching_cubes(SDF2, GV2, grid_resolution, grid_resolution, grid_resolution, 0, V3, F3);
//
//
//    // --- 3. 准备布尔运算 (翻转法线) ---
//   // Marching Cubes生成的法线通常朝内，我们需要翻转它
////    F.col(1).swap(F.col(2));
////    std::cout << "Normals of Model B flipped." << std::endl;
////
////    // --- 4. 执行布尔求差 (A - B) ---
////    Eigen::MatrixXd VC;
////    Eigen::MatrixXi FC;
////    igl::copyleft::cgal::mesh_boolean(VA, FA, V, F, igl::MESH_BOOLEAN_TYPE_UNION, VC, FC);
//////    igl::copyleft::cgal::mesh_boolean(VA, FA, V, F, igl::MESH_BOOLEAN_TYPE_MINUS, VC, FC);
////    std::cout << "Boolean difference (A - B) operation completed." << std::endl;
////
////    // --- 5. 保存结果 ---
////    // 保存高斯球模型
////    igl::write_triangle_mesh("gaussian_spheres_model.obj", V, F);
////    std::cout << "Model B saved to 'gaussian_spheres_model.obj'" << std::endl;
////
////    // 保存布尔求差后的模型
////    igl::write_triangle_mesh("boolean_difference_result.obj", VC, FC);
////    std::cout << "Boolean result saved to 'boolean_difference_result.obj'" << std::endl;
//
//
//    // =========================================================================
//    // 5. 使用 libigl Viewer 显示结果
//    // =========================================================================
//    std::cout << "启动 libigl 查看器..." << std::endl;
//    igl::opengl::glfw::Viewer viewer;
//
//    //Eigen::MatrixXi F_flip = F;
//    //F_flip.col(0).swap(F_flip.col(1));  // 交换0和1列
//    //F = F_flip;
//    // 设置主网格 (高斯球)
//    viewer.data().set_mesh(VA, FA);
//    viewer.data().set_face_based(true); // 使用基于面的法线，获得更平滑的外观
//    viewer.data().show_lines = false;   // 不显示网格线
//   //viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色
//
//    // 添加辅助点 (高斯核的中心)，设置为红色
//   // viewer.data().add_points(kernel_points, Eigen::RowVector3d(1, 0, 0));
//    viewer.data().point_size = 10; // 让点更显眼
//
//    viewer.launch();
//
//
//
//    //std::cout << "启动 libigl 查看器..." << std::endl;
//    //igl::opengl::glfw::Viewer viewer2;
//
//    ////Eigen::MatrixXi F_flip = F;
//    ////F_flip.col(0).swap(F_flip.col(1));  // 交换0和1列
//    ////F = F_flip;
//    //// 设置主网格 (高斯球)
//    //viewer2.data().set_mesh(V3, F3);
//    //viewer2.data().set_face_based(true); // 使用基于面的法线，获得更平滑的外观
//    //viewer2.data().show_lines = false;   // 不显示网格线
//    ////viewer.data().set_colors(Eigen::RowVector3d(0.8, 0.7, 0.2)); // 设置一个漂亮的蓝色
//
//    // // 添加辅助点 (高斯核的中心)，设置为红色
//    //// viewer.data().add_points(kernel_points, Eigen::RowVector3d(1, 0, 0));
//    //viewer2.data().point_size = 10; // 让点更显眼
//
//    //viewer2.launch();
//
//
//    return 0;
//}

#include <igl/marching_cubes.h>
#include <igl/signed_distance.h>
#include <igl/read_triangle_mesh.h>
#include <igl/voxel_grid.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>


int main(int argc, char* argv[])
{
    using namespace Eigen;
    using namespace std;
    using namespace igl;
    MatrixXi F;
    MatrixXd V;
    // Read in inputs as double precision floating point meshes
    read_triangle_mesh(
        TUTORIAL_SHARED_PATH "/armadillo.obj", V, F);
    cout << "Creating grid..." << endl;
    // number of vertices on the largest side
    const int s = 100;
    // create grid
    MatrixXd GV;
    Eigen::RowVector3i res;
    igl::voxel_grid(V, 0, s, 1, GV, res);

    // compute values
    cout << "Computing distances..." << endl;
    VectorXd S, B;
    {
        VectorXi I;
        MatrixXd C, N;
        signed_distance(GV, V, F, SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C, N);
        // Convert distances to binary inside-outside data --> aliasing artifacts
        B = S;
        for_each(B.data(), B.data() + B.size(), [](double& b) {b = (b > 0 ? 1 : (b < 0 ? -1 : 0)); });
    }
    cout << "Marching cubes..." << endl;
    MatrixXd SV, BV;
    MatrixXi SF, BF;
    igl::marching_cubes(S, GV, res(0), res(1), res(2), 0, SV, SF);
    igl::marching_cubes(B, GV, res(0), res(1), res(2), 0, BV, BF);

    cout << R"(Usage:
'1'  Show original mesh.
'2'  Show marching cubes contour of signed distance.
'3'  Show marching cubes contour of indicator function.
)";
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(SV, SF);
    viewer.callback_key_down =
        [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod)->bool
        {
            switch (key)
            {
            default:
                return false;
            case '1':
                viewer.data().clear();
                viewer.data().set_mesh(V, F);
                break;
            case '2':
                viewer.data().clear();
                viewer.data().set_mesh(SV, SF);
                break;
            case '3':
                viewer.data().clear();
                viewer.data().set_mesh(BV, BF);
                break;
            }
            viewer.data().set_face_based(true);
            return true;
        };
    viewer.launch();
}
