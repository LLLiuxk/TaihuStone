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
    std::cout << "--------------------------" << std::endl;


    std::string input_path = "D:/VSprojects/TaihuStone/model/" + input_file + ".stl";

    ModelGenerator mg(input_path);
    //mg.test_item();
    mg.model_porous_structure();
    if (figure_show)
        mg.show_model();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Model generation completed in " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
}
