#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cstdint>
#include <utility>

// ============================================================================
// TPMS Pipeline: from TPMS implicit field to skeleton graph + VP analysis
// ============================================================================
// Pipeline steps:
//   1. Generate regular / perturbed TPMS implicit field
//   2. Discretize into 3D voxel grid
//   3. Extract largest connected void region
//   4. 3D skeleton extraction (topology-preserving thinning)
//   5. Prune short branches
//   6. skeleton -> graph (junction clustering + node classification)
//   7. Compute VP (Visual Permeability) on graph
//   8. Export visualization data
// ============================================================================

// ---------- type aliases ----------
using Vector3d = Eigen::Vector3d;
using Vector3i = Eigen::Vector3i;
using MatrixXd = Eigen::MatrixXd;
using MatrixXi = Eigen::MatrixXi;

// ---------- TPMS parameters ----------
struct TpmsParams
{
    // spatial frequency perturbation: alpha_* near 1.0, with small random perturbation
    double alpha_x = 1.0;
    double alpha_y = 1.0;
    double alpha_z = 1.0;

    // phase offsets
    double phase_x = 0.0;
    double phase_y = 0.0;
    double phase_z = 0.0;

    // voxel resolution (per axis)
    int    resolution = 64;

    // physical bounding box: [-half_bbox, +half_bbox]^3
    double half_bbox = 3.141592653589793;   // ~pi: one full Gyroid period

    // iso-surface value (f=0 is the classic Gyroid iso-surface)
    double iso_value = 0.0;

    // whether this is a perturbed version
    bool is_perturbed = false;

    // voxel spacing = 2*half_bbox / (resolution-1)
    double spacing() const { return 2.0 * half_bbox / double(resolution - 1); }

    // total grid points
    int grid_total() const { return resolution * resolution * resolution; }
};

// ---------- 3D integer voxel field ----------
struct VoxelField3D
{
    int nx = 0, ny = 0, nz = 0;
    std::vector<uint8_t> data;   // 0=empty, 1=filled (belongs to void region)

    int total() const { return nx * ny * nz; }

    inline int idx(int x, int y, int z) const
    {
        return x + nx * (y + ny * z);
    }

    inline bool in_bounds(int x, int y, int z) const
    {
        return x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz;
    }

    inline uint8_t at(int x, int y, int z) const { return data[idx(x, y, z)]; }
    inline uint8_t& at(int x, int y, int z) { return data[idx(x, y, z)]; }
    inline uint8_t at_idx(int i) const { return data[i]; }
    inline uint8_t& at_idx(int i) { return data[i]; }

    inline Vector3d world_pos(int x, int y, int z,
        const Vector3d& origin, double spacing) const
    {
        return Vector3d(
            origin.x() + x * spacing,
            origin.y() + y * spacing,
            origin.z() + z * spacing
        );
    }

    inline void coord(int flat_idx, int& x, int& y, int& z) const
    {
        x = flat_idx % nx;
        int rem = flat_idx / nx;
        y = rem % ny;
        z = rem / ny;
    }
};

// ---------- skeleton voxel (with distance field for classification) ----------
struct SkeletonVoxel
{
    int  idx_flat = -1;          // flat index in voxel field
    int  x = 0, y = 0, z = 0;
    int  degree = 0;             // number of skeleton neighbors in 26-neighborhood
    double dist_to_surface = 0;  // distance to nearest TPMS wall (voxel units)
    Vector3d world_pos;
};

// ---------- graph node type ----------
enum class GraphNodeType : uint8_t
{
    JUNCTION = 0,   // branch junction (clustered skeleton voxels with degree >= 3)
    ENDPOINT = 1,   // endpoint (degree == 1)
    VIRTUAL = 2    // virtual node inserted in long edges
};

// ---------- graph node location classification ----------
enum class NodeLocation : uint8_t
{
    INTERIOR = 0,   // deep in skeleton, far from TPMS wall
    SURFACE = 1    // near TPMS wall (< dist_threshold)
};

// ---------- graph node ----------
struct GraphNode
{
    int         id = -1;
    GraphNodeType type = GraphNodeType::JUNCTION;
    NodeLocation  location = NodeLocation::INTERIOR;
    Vector3d    position = Vector3d::Zero();
    double      dist_to_surface = 0.0;
    std::vector<int> neighbors;   // indices of adjacent graph edges
};

// ---------- graph edge ----------
struct GraphEdge
{
    int  id = -1;
    int  node_a = -1;
    int  node_b = -1;
    double length_world = 0.0;               // world-space length
    std::vector<Vector3d> polyline;           // skeleton polyline (ordered, includes endpoints)
};

// ---------- VP weights ----------
struct VpWeights
{
    double w_angle = 0.25;
    double w_length = 0.25;
    double w_location = 0.25;
    double w_direction = 0.25;   // weight for S_horiz
};

// ---------- VP computation result ----------
struct VpResult
{
    double total_vp = 0.0;               // average VP over all nodes
    std::vector<double> node_vp;          // VP per node
    std::vector<std::vector<int>> node_best_path; // best path per node
};

// ---------- pipeline statistics ----------
struct PipelineStats
{
    // TPMS
    int    resolution = 0;
    double bbox_half = 0.0;

    // Void region
    int    total_void_voxels = 0;
    int    largest_void_voxels = 0;

    // Skeleton
    int    skeleton_voxels_initial = 0;
    int    skeleton_voxels_pruned = 0;

    // Graph
    int    graph_node_count = 0;
    int    graph_edge_count = 0;
    int    surface_node_count = 0;
    int    interior_node_count = 0;
    double avg_degree = 0.0;

    // VP
    double final_vp = 0.0;

    // Meta
    bool   is_perturbed = false;
};

// ============================================================================
// Main Pipeline class
// ============================================================================
class TpmsPipeline
{
public:
    TpmsPipeline() = default;

    // ---------- parameter setters ----------
    void setParams(const TpmsParams& p) { params_ = p; }
    void setVpWeights(const VpWeights& w) { vp_weights_ = w; }

    // ---------- run full pipeline ----------
    // steps: bitmask, 0=all. Can be used for step-by-step debugging.
    PipelineStats run(uint32_t steps_mask = 0);

    // ---------- individual step interfaces (for debugging/replacement) ----------
    void step1_generateTpmsField();
    void step2_voxelize();
    void step3_extractLargestVoidRegion();
    void step4_skeletonize();
    void step5_pruneSkeleton(double min_branch_length_voxels = 3.0);
    void step6_buildGraph(double junction_cluster_radius = 1.8, double max_edge_no_node = 8.0);
    void step7_computeVP();
    void step8_exportData(const std::string& output_prefix);

    // ---------- query interfaces ----------
    const PipelineStats& stats()      const { return stats_; }
    const std::vector<GraphNode>& graphNodes() const { return graph_nodes_; }
    const std::vector<GraphEdge>& graphEdges() const { return graph_edges_; }
    const VpResult& vpResult()   const { return vp_result_; }
    const std::vector<double>& tpmsField() const { return tpms_field_; }

private:
    // ---- data ----
    TpmsParams  params_;
    VpWeights   vp_weights_;

    // TPMS continuous field (flattened, size = res^3)
    std::vector<double> tpms_field_;

    // voxelized void mask
    VoxelField3D void_mask_;         // 1 = void (f < iso - margin)
    VoxelField3D largest_void_;      // largest connected component

    // skeleton
    std::vector<uint8_t> skeleton_data_;   // flattened, 1 = skeleton voxel
    std::vector<SkeletonVoxel> skeleton_voxels_;

    // graph
    std::vector<GraphNode> graph_nodes_;
    std::vector<GraphEdge> graph_edges_;

    // VP
    VpResult vp_result_;

    // statistics
    PipelineStats stats_;

    // ---- internal helpers ----
    Vector3d grid_pos(int x, int y, int z) const;

    // 6-neighbor & 26-neighbor offsets
    static const int k6[6][3];
    static const int k26[26][3];

    // connected components (6-neighbor for void)
    std::vector<std::vector<int>> extract_components_6(const VoxelField3D& field);

    // 3D thinning helpers
    int  count_neighbors_26(int flat_idx, const VoxelField3D& field) const;
    int  count_neighbors_6(int flat_idx, const VoxelField3D& field) const;
    bool is_simple_point(int flat_idx, const VoxelField3D& field) const;
    bool is_endpoint(int flat_idx) const;
    void distance_transform(const VoxelField3D& src, std::vector<double>& dist) const;

    // connected components (26-neighbor, for skeleton)
    std::vector<std::vector<int>> extract_components_26(const std::vector<uint8_t>& data,
        int nx, int ny, int nz) const;

    // graph construction helper
    void trace_skeleton_path(int start_flat, int start_x, int start_y, int start_z,
        const std::vector<int>& node_of_flat,
        std::vector<Vector3d>& path, int& reached_node) const;

    // VP helpers
    double compute_path_vp(const std::vector<Vector3d>& path_points,
        const std::vector<bool>& is_surface) const;
    Eigen::Vector3d computePrincipalDirection(const std::vector<Eigen::Vector3d>& pts) const;
    double abs_angle_deg(const Vector3d& v1, const Vector3d& v2) const;
};

// ============================================================================
// Convenience function: generate two TPMS cases and run full pipeline
// ============================================================================
void runTpmsComparison(int resolution, const std::string& output_dir);
