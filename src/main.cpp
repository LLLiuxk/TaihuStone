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

 //   Eigen::Vector3d z(0, 0, 1);
 //   /*std::vector<Eigen::Vector3d> path_points = {
 //       Eigen::Vector3d(0.153543,  0.011811, -0.279528),
 //       Eigen::Vector3d(0.0748031, -0.0984252, -0.248031),
 //       Eigen::Vector3d(0.0275591, -0.0433071, -0.161417),
 //       Eigen::Vector3d(-0.0748031,  0.0354331, -0.129921),
 //       Eigen::Vector3d(-0.0748031, -0.0433071, -0.216535),
 //       Eigen::Vector3d(-0.114173,  -0.208661, -0.279528)
	//};*/
 //   /*std::vector<Eigen::Vector3d> path_points = {
 //       Eigen::Vector3d(0, 0, 0),
 //       Eigen::Vector3d(3, 0, -2),
 //       Eigen::Vector3d(5, 0, 2),
 //       Eigen::Vector3d(7, 0, -2),
 //       Eigen::Vector3d(9, 0, 2),
 //       Eigen::Vector3d(10, 0, 0),
 //       Eigen::Vector3d(1, 0, 2)
 //   };*/
 //   std::vector<Eigen::Vector3d> path_points = {
 //       Eigen::Vector3d(0, 0, 0),
 //       Eigen::Vector3d(3, 3, 1),
 //       Eigen::Vector3d(5, 1, 3),
 //       Eigen::Vector3d(3, 3, 1),
 //       Eigen::Vector3d(5, 1, 3),
 //       Eigen::Vector3d(3, 3, 1),
 //       Eigen::Vector3d(5, 1, 3),
 //       Eigen::Vector3d(3, 3, 1),
 //       Eigen::Vector3d(5, 1, 3),
 //       Eigen::Vector3d(7, -2, 0),
 //       Eigen::Vector3d(10, 0, 0)
 //   };
 //   Eigen::Vector3d dir = computePrincipalDirection(path_points);
	//cout << "Principal Direction: " << dir.transpose() << endl;
 //   double S_horiz = 1.0 - std::abs(dir.dot(z));
	//cout << "S_horiz: " << S_horiz << endl;


    std::string configFileName = "default"; // 用于日志显示
    if (argc > 1) {
        configFileName = argv[1];
        std::cout << "Loading configuration from: " << configFileName << std::endl;
        if (!loadParameters(configFileName)) {
            std::cerr << "Failed to load parameters, exiting." << std::endl;
            return -1;
        }
    }
    else {
        std::cout << "No configuration file provided. Using default hardcoded parameters." << std::endl;
    }
    std::cout << "--- Current Parameters ---" << std::endl;
    std::cout << "Input File: " << input_file << std::endl;
    std::cout << "PoresNum: " << PoresNum << std::endl;
    std::cout << "compare_show: " << compare_show << std::endl;
    std::cout << "topo_optimize: " << topo_optimize << std::endl;
    std::cout << "figure_show: " << figure_show << std::endl;
    std::cout << "Trans_thres: " << Trans_thres << std::endl;
    std::cout << "--------------------------" << std::endl;


    std::string input_path = "D:/VSprojects/TaihuStone/model/" + input_file + ".stl";

    ModelGenerator mg(input_path);
    //mg.test_item();
    mg.model_porous_structure();
    if (scr_figure)
        mg.show_model();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Model generation completed in " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
}
