#include <igl/marching_cubes.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/cgal/mesh_boolean.h>

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <random>
#include <cstdio>
#include <array>


#include "Tool.h"
#include "modelGen.h"
#include "TpmsPipeline.h"

static void graphToTubeMesh(
    const std::vector<GraphNode>& nodes,
    const std::vector<GraphEdge>& edges,
    Eigen::MatrixXd& V_out, Eigen::MatrixXi& F_out,
    double tube_radius, int n_sides)
{
    std::vector<Eigen::Vector3d> V;
    std::vector<Eigen::Vector3i> F;

    const double PI = 3.141592653589793;

    for (const auto& e : edges) {
        const auto& poly = e.polyline;
        for (size_t i = 0; i + 1 < poly.size(); ++i) {
            Eigen::Vector3d p0 = poly[i];
            Eigen::Vector3d p1 = poly[i + 1];
            Eigen::Vector3d dir = (p1 - p0).normalized();

            Eigen::Vector3d ref = (std::abs(dir.z()) < 0.999)
                ? Eigen::Vector3d(0, 0, 1) : Eigen::Vector3d(1, 0, 0);
            Eigen::Vector3d u = dir.cross(ref).normalized();
            Eigen::Vector3d v = dir.cross(u).normalized();

            int base = (int)V.size();
            for (int s = 0; s < n_sides; ++s) {
                double ang = 2.0 * PI * s / n_sides;
                Eigen::Vector3d offset = tube_radius * (std::cos(ang) * u + std::sin(ang) * v);
                V.push_back(p0 + offset);
            }
            for (int s = 0; s < n_sides; ++s) {
                double ang = 2.0 * PI * s / n_sides;
                Eigen::Vector3d offset = tube_radius * (std::cos(ang) * u + std::sin(ang) * v);
                V.push_back(p1 + offset);
            }

            for (int s = 0; s < n_sides; ++s) {
                int s1 = (s + 1) % n_sides;
                int a = base + s;
                int b = base + s1;
                int c = base + n_sides + s1;
                int d = base + n_sides + s;
                F.push_back(Eigen::Vector3i(a, b, c));
                F.push_back(Eigen::Vector3i(a, c, d));
            }
        }
    }

    for (const auto& nd : nodes) {
        double r = tube_radius * 1.5;
        Eigen::Vector3d c = nd.position;
        int base = (int)V.size();
        V.push_back(c + Eigen::Vector3d(r, 0, 0));
        V.push_back(c + Eigen::Vector3d(-r, 0, 0));
        V.push_back(c + Eigen::Vector3d(0, r, 0));
        V.push_back(c + Eigen::Vector3d(0, -r, 0));
        V.push_back(c + Eigen::Vector3d(0, 0, r));
        V.push_back(c + Eigen::Vector3d(0, 0, -r));
        int b = base;
        F.push_back(Eigen::Vector3i(b + 0, b + 2, b + 4));
        F.push_back(Eigen::Vector3i(b + 2, b + 1, b + 4));
        F.push_back(Eigen::Vector3i(b + 1, b + 3, b + 4));
        F.push_back(Eigen::Vector3i(b + 3, b + 0, b + 4));
        F.push_back(Eigen::Vector3i(b + 0, b + 5, b + 2));
        F.push_back(Eigen::Vector3i(b + 2, b + 5, b + 1));
        F.push_back(Eigen::Vector3i(b + 1, b + 5, b + 3));
        F.push_back(Eigen::Vector3i(b + 3, b + 5, b + 0));
    }

    V_out.resize(V.size(), 3);
    for (size_t i = 0; i < V.size(); ++i) V_out.row(i) = V[i];
    F_out.resize(F.size(), 3);
    for (size_t i = 0; i < F.size(); ++i) F_out.row(i) = F[i];
}

int main(int argc, char* argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();
    int chosen_fun = 1;
    // ========================================================================
    // TPMS + Skeleton preview
    //   1) TPMS field → marching cubes → STL (surface)
    //   2) Voxelize → skeletonize → build graph → tube mesh STL (skeleton)
    // ========================================================================
    if (chosen_fun == 0)
    {
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
        std::cout << "Max_edge_limit: " << Max_edge_limit << std::endl;
        std::cout << "compare_show: " << compare_show << std::endl;
        std::cout << "topo_optimize: " << topo_optimize << std::endl;
        std::cout << "figure_show: " << figure_show << std::endl;
        std::cout << "Trans_thres: " << Trans_thres << std::endl;
        std::cout << "--------------------------" << std::endl;


        std::string input_path = "D:/VSprojects/TaihuStone/model/" + input_file + ".stl";

        ModelGenerator mg(input_path);
        //mg.test_item(); //smooth
        mg.model_porous_structure();
        if (scr_figure)
            mg.show_model();

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Model generation completed in " << duration.count() / 1000.0 << " s" << std::endl;
    }
    else if (chosen_fun == 1)
    {
        std::string out_dir = "D:/VSprojects/TaihuStone/tpms_output";
        int res = 64;

        double sp = 2.0 * M_PI / (res - 1);
        double tube_r = sp * 0.35;

        // ---- Case 1: Regular Gyroid ----
        std::cout << "\n========== Regular Gyroid (alpha=2.5) ==========\n";
        {
            TpmsParams p;
            p.resolution = res;
            p.alpha_x = 2.5; p.alpha_y = 2.5; p.alpha_z = 2.5;
            p.phase_x = 0.0; p.phase_y = 0.0; p.phase_z = 0.0;
            p.is_perturbed = false;

            TpmsPipeline pipe;
            pipe.setParams(p);
            PipelineStats s = pipe.run(63);

            int R = p.resolution;
            double ox = -p.half_bbox;
            const auto& field = pipe.tpmsField();

            Eigen::MatrixXd GV(R * R * R, 3);
            Eigen::VectorXd S(R * R * R);
            for (int z = 0; z < R; ++z)
                for (int y = 0; y < R; ++y)
                    for (int x = 0; x < R; ++x) {
                        int idx = x + R * (y + R * z);
                        GV.row(idx) = Eigen::Vector3d(ox + x * sp, ox + y * sp, ox + z * sp);
                        S(idx) = field[idx];
                    }

            Eigen::MatrixXd V_surf, V_skel;
            Eigen::MatrixXi F_surf, F_skel;

            MarchingCubes(S, GV, R, R, R, 0.0, V_surf, F_surf);
            saveMesh(out_dir + "/tpms_regular_preview.stl", V_surf, F_surf);

            graphToTubeMesh(pipe.graphNodes(), pipe.graphEdges(), V_skel, F_skel, tube_r, 8);
            saveMesh(out_dir + "/tpms_regular_skeleton.stl", V_skel, F_skel);

            std::cout << "  Surface: " << V_surf.rows() << " verts, " << F_surf.rows() << " faces\n";
            std::cout << "  Skeleton: " << s.skeleton_voxels_pruned << " voxels -> "
                << pipe.graphNodes().size() << " nodes, " << pipe.graphEdges().size() << " edges"
                << " (tube mesh: " << V_skel.rows() << " verts, " << F_skel.rows() << " faces)\n";
            std::cout << "  Avg degree: " << s.avg_degree << "\n";
        }

        // ---- Case 2: Perturbed Gyroid ----
        std::cout << "\n========== Perturbed Gyroid (alpha ~2.5-3.5) ==========\n";
        {
            std::random_device rd;
            std::mt19937 rng(rd());
            std::uniform_real_distribution<double> alpha_noise(2.0, 3.0);
            std::uniform_real_distribution<double> phase_noise(-0.25, 0.25);

            TpmsParams p;
            p.resolution = res;
            p.alpha_x = alpha_noise(rng); p.alpha_y = alpha_noise(rng); p.alpha_z = alpha_noise(rng);
            p.phase_x = phase_noise(rng); p.phase_y = phase_noise(rng); p.phase_z = phase_noise(rng);
            p.is_perturbed = true;
            std::cout << "  alpha=(" << p.alpha_x << "," << p.alpha_y << "," << p.alpha_z << ")\n";
            std::cout << "  phase=(" << p.phase_x << "," << p.phase_y << "," << p.phase_z << ")\n";

            TpmsPipeline pipe;
            pipe.setParams(p);
            PipelineStats s = pipe.run(63);

            int R = p.resolution;
            double ox = -p.half_bbox;
            const auto& field = pipe.tpmsField();

            Eigen::MatrixXd GV(R * R * R, 3);
            Eigen::VectorXd S(R * R * R);
            for (int z = 0; z < R; ++z)
                for (int y = 0; y < R; ++y)
                    for (int x = 0; x < R; ++x) {
                        int idx = x + R * (y + R * z);
                        GV.row(idx) = Eigen::Vector3d(ox + x * sp, ox + y * sp, ox + z * sp);
                        S(idx) = field[idx];
                    }

            Eigen::MatrixXd V_surf, V_skel;
            Eigen::MatrixXi F_surf, F_skel;

            MarchingCubes(S, GV, R, R, R, 0.0, V_surf, F_surf);
            saveMesh(out_dir + "/tpms_perturbed_preview.stl", V_surf, F_surf);

            graphToTubeMesh(pipe.graphNodes(), pipe.graphEdges(), V_skel, F_skel, tube_r, 8);
            saveMesh(out_dir + "/tpms_perturbed_skeleton.stl", V_skel, F_skel);

            std::cout << "  Surface: " << V_surf.rows() << " verts, " << F_surf.rows() << " faces\n";
            std::cout << "  Skeleton: " << s.skeleton_voxels_pruned << " voxels -> "
                << pipe.graphNodes().size() << " nodes, " << pipe.graphEdges().size() << " edges"
                << " (tube mesh: " << V_skel.rows() << " verts, " << F_skel.rows() << " faces)\n";
            std::cout << "  Avg degree: " << s.avg_degree << "\n";
        }

        std::cout << "\n=== Preview done. Files saved to " << out_dir << " ===\n";
        std::cout << "    *_preview.stl = TPMS surface\n";
        std::cout << "    *_skeleton.stl = skeleton graph (tube mesh, rotatable)\n";
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Completed in " << duration.count() / 1000.0 << " s" << std::endl;

    return 0;
}
