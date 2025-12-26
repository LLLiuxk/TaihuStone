#pragma once
#include <igl/readOBJ.h>
#include <igl/signed_distance.h>

#include <igl/marching_cubes.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/voxel_grid.h>
#include <igl/opengl/glfw/Viewer.h>

#include <filesystem>
#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <queue>
#include <random>
#include "FastNoiseLite.h"
#include "globalPara.h" 
#include "selfSupVoxel.h"

#define M_PI 3.1415926
#define INF std::numeric_limits<double>::infinity()
using namespace std;
using namespace Eigen;


static FastNoiseLite g_kernel_noise;
static FastNoiseLite g_field_noise;

//static void init_noise()
//{
//    static bool initialized = false;
//    if (initialized) return;
//
//    g_kernel_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
//    g_kernel_noise.SetFractalType(FastNoiseLite::FractalType_FBm);
//    g_kernel_noise.SetFractalOctaves(3);        // 低频
//    g_kernel_noise.SetFractalLacunarity(2.0f);
//    g_kernel_noise.SetFractalGain(0.5f);
//    g_kernel_noise.SetFrequency(0.8f);          // 核心参数：空间尺度
//    initialized = true;
//}

static void init_field_noise()
{
    static bool initialized = false;
    if (initialized) return;

    g_field_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    g_field_noise.SetFractalType(FastNoiseLite::FractalType_FBm);
    g_field_noise.SetFractalOctaves(4);  //高频
    g_field_noise.SetFractalGain(0.55f);
    g_field_noise.SetFrequency(1.0f); // 空间尺度由外部控制

    initialized = true;
}

double Mesh2SDF(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& GV, Eigen::VectorXd& S, Eigen::Vector3d& bb_min, Eigen::Vector3d& bb_max);
bool saveMesh(std::string filename, Eigen::MatrixXd V, Eigen::MatrixXi F);

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
double distance(Vector3d v1, Vector3d v2);
double squared_distance(Vector3d v1, Vector3d v2);

//Bernstein基函数 
double bernstein_basis(int i, int n, double t);

//load files
bool exportSDF(Eigen::VectorXd& sdf, std::string& filename);

bool align_models_with_pca(const std::string& model1_path, const std::string& model2_path, const std::string& output_path);

//show_result 
void show_path(std::vector<int> path);

//kinds of check
void geometry_analyzer(Eigen::VectorXd SDF, int resolution, double thres_degree, int& overhang_count, int& floating_count, std::vector<uint8_t>& overhang_mask, std::vector<uint8_t>& floating_mask);
Vector3d computeGradient(int x, int y, int z, int res, Eigen::VectorXd SDF);
void getCoord(int idx, int res, int& x, int& y, int& z);

// pca point cloud
Eigen::Vector3d computePrincipalDirection(const std::vector<Eigen::Vector3d>& points);

//TO
double smoothHeaviside(double s, double eps = 0.1);  // Heaviside 平滑宽度（建议 = 1~2 个体素尺寸）
double hardTrans(double s, double iso);

VoxelGrid SDFtoVoxel(Eigen::VectorXd& sdf, Eigen::Vector3d minBox, Eigen::Vector3d maxBox, int nx, int ny, int nz); 

void saveVoxelToRaw(std::string filename, VoxelGrid& grid);

void saveVoxelToVTK(std::string filename, VoxelGrid& grid);

void saveVoxelGridAsNPY(std::vector<double>& voxel_grid, int res, std::string& filename);

void writeAdjacencyListToFile(const std::vector<std::vector<int>>& adjList, const std::string& filename);

std::vector<std::vector<int>> readAdjacencyListFromFile(const std::string& filename);

//string
std::string to_string_pre(double value, int precision = 1);

//sort
template <typename T1, typename T2>
void sort_min2max(std::vector<std::pair<T1, T2>>& vec)
{
    std::stable_sort(vec.begin(), vec.end(),
        [](const std::pair<T1, T2>& a, const std::pair<T1, T2>& b)
        {
            return a.second < b.second;
        });
}

//double add_iso_surface_noise(
//    const Eigen::Vector3d& p,
//    double sdf_value,
//    double band_width,
//    double noise_amplitude,
//    double spatial_frequency);

void add_noise_near_isosurface(
    Eigen::VectorXd& S,              // 标量场（in-place 修改）
    const Eigen::MatrixXd& GV,        // 体素坐标
    double iso_value,                 // MC 等值面
    double noise_amplitude,            // 噪声幅值（建议 0.05 ~ 0.3 * 场尺度）
    double band_width,                 // 等值面带宽
    double spatial_frequency           // 噪声空间频率
);

//cal sigma
double sigma_value(double v, double n, double w, double iso);
