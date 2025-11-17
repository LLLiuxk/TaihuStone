#pragma once
#include <igl/readOBJ.h>
#include <igl/signed_distance.h>

#include <igl/marching_cubes.h>
#include <igl/read_triangle_mesh.h>
#include <igl/voxel_grid.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/write_triangle_mesh.h>

#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

#include "globalPara.h" 

#define M_PI 3.1415926
using namespace std;
using namespace Eigen;

struct ErrorPoint {
    int x, y, z;
    double value1, value2, error;

};
void Mesh2SDF(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& GV, Eigen::VectorXd& S);


// SDF平滑并集。k越大，平滑效果越小，趋近于普通并集
double smooth_UnionSDF(double sdf1, double sdf2, double k);

// SDF平滑交集
double smooth_IntersecSDF(double sdf1, double sdf2, double k);

// SDF布尔运算：并集
double unionSDF(double sdf1, double sdf2);

// SDF布尔运算：交集
double intersectionSDF(double sdf1, double sdf2);
// SDF布尔运算：差集 (A - B)
double differenceSDF(double sdf1, double sdf2);

void MarchingCubes(Eigen::VectorXd& S, Eigen::MatrixXd& GV, int nx, int ny, int nz, double isovalue, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

void view_model(Eigen::MatrixXd V1, Eigen::MatrixXi F1);
void view_two_models(Eigen::MatrixXd V1, Eigen::MatrixXi F1, Eigen::MatrixXd V2, Eigen::MatrixXi F2, Eigen::RowVector3d shift = RowVector3d(0.0, 0.0, 0.0));
void view_three_models(Eigen::MatrixXd V1, Eigen::MatrixXi F1, Eigen::MatrixXd V2, Eigen::MatrixXi F2, Eigen::MatrixXd V3, Eigen::MatrixXi F3, Eigen::RowVector3d shift = RowVector3d(0.0, 0.0, 0.0));

int  single_component(Eigen::MatrixXd V, Eigen::MatrixXi F);

//三角运算
double abs_angle(Vector3d v1, Vector3d v2);

//Bernstein基函数 
double bernstein_basis(int i, int n, double t);

//load files
bool exportSDF(Eigen::VectorXd& sdf, std::string& filename);

bool align_models_with_pca(const std::string& model1_path, const std::string& model2_path, const std::string& output_path);
 
