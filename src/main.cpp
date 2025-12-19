#include <igl/marching_cubes.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/cgal/mesh_boolean.h>

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <random>
#include <cstdio>

#include "Tool.h"
#include "modelGen.h"



int main(int argc, char* argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();

    std::string input_file = "D:/VSprojects/TaihuStone/model/RockSet03r.stl";
    //std::string input_file = "D:/VSprojects/TaihuStone/model/test_case.stl";

    ModelGenerator mg(input_file);
   
    //mg.show_model();


    //----------test---------------
    //std::string filename = "result/36-14.txt";
    //std::vector<std::vector<int>> adj = readAdjacencyListFromFile(filename);
    //PathQuery p_bfs(PoresNum, adj, 36);
    //int s1 = 14;
    //int s2 = 6;

    //std::vector<int> path1 = p_bfs.query_path(s1);
    //std::vector<int> path2 = p_bfs.query_path(s2);
    //cout << path1.size() << "   " << path2.size() << endl;


 //   std::string input_file_ = "D:/VSprojects/TaihuStone/model/gaussian_pores150.stl";
 //   Eigen::MatrixXd V_ini; //初始网格顶点
 //   Eigen::MatrixXi F_ini; // 初始网格面片
 //   Eigen::VectorXd SDF_ini;           // 原始/已处理的SDF网格值
 //   Eigen::MatrixXd GV;            // 网格点坐标
 //   if (!igl::read_triangle_mesh(input_file_, V_ini, F_ini)) {
 //       std::cerr << "Error: Could not load model A." << std::endl;
 //   }
 //   Mesh2SDF(V_ini, F_ini, GV, SDF_ini);
 //   Eigen::MatrixXd V_t; //输出网格顶点
 //   Eigen::MatrixXi F_t; // 输出网格面片
 //   MarchingCubes(SDF_ini, GV, Resolution, Resolution, Resolution, 0, V_t, F_t);  //gaussian combined with tubes
 //   view_model(V_t, F_t);

	////VoxelGrid grid_ = SDFtoVoxel(SDF_ini, GV.row(0), GV.row(GV.rows() - 1), Resolution, Resolution, Resolution, 1.0);
	////saveVoxelToRaw( "D:/VSprojects/TaihuStone/model/raw/voxelized_model.raw", grid_);
	////saveVoxelToVTK("D:/VSprojects/TaihuStone/model/vtk/voxelized_model.vtk", grid_);
 //   //  保存体素化结果为NPY格式
 //   std::cout << "开始保存体素化结果..." << std::endl;
 //   // 创建二值体素网格 (1=实心, 0=空心)
 //   std::vector<uint8_t> voxel_grid(Resolution * Resolution * Resolution, 0);

 //   for (int i = 0; i < SDF_ini.size(); ++i) {
 //       voxel_grid[i] = (SDF_ini(i) < 0) ? 1 : 0;
 //   }

 //   // 保存为NPY格式
 //   std::string npy_filename = "D:/VSprojects/TaihuStone/model/npy/voxelized_model.npy";
 //   saveVoxelGridAsNPY(voxel_grid, Resolution, npy_filename);

	// ------------------------------------- view two models after alignment -------------------------------------
 
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Model generation completed in " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
}
