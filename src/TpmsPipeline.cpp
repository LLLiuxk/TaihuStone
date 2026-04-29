#include "TpmsPipeline.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <stack>
#include <set>
#include <map>
#include <chrono>
#include <cassert>

using namespace std;

// ============================================================================
// 静态常量: 邻域偏移
// ============================================================================
const int TpmsPipeline::k6[6][3] = {
    {+1,0,0},{-1,0,0},{0,+1,0},{0,-1,0},{0,0,+1},{0,0,-1}
};

const int TpmsPipeline::k26[26][3] = {
    {-1,-1,-1},{-1,-1, 0},{-1,-1, 1},
    {-1, 0,-1},{-1, 0, 0},{-1, 0, 1},
    {-1, 1,-1},{-1, 1, 0},{-1, 1, 1},
    { 0,-1,-1},{ 0,-1, 0},{ 0,-1, 1},
    { 0, 0,-1},            { 0, 0, 1},
    { 0, 1,-1},{ 0, 1, 0},{ 0, 1, 1},
    { 1,-1,-1},{ 1,-1, 0},{ 1,-1, 1},
    { 1, 0,-1},{ 1, 0, 0},{ 1, 0, 1},
    { 1, 1,-1},{ 1, 1, 0},{ 1, 1, 1}
};

// ============================================================================
// 辅助函数
// ============================================================================
Vector3d TpmsPipeline::grid_pos(int x, int y, int z) const {
    double sp = params_.spacing();
    double ox = -params_.half_bbox;
    return Vector3d(ox + x * sp, ox + y * sp, ox + z * sp);
}

double TpmsPipeline::abs_angle_deg(const Vector3d& v1, const Vector3d& v2) const {
    double dot = v1.normalized().dot(v2.normalized());
    dot = clamp(dot, -1.0, 1.0);
    return acos(dot) * 180.0 / 3.141592653589793;
}

Eigen::Vector3d TpmsPipeline::computePrincipalDirection(const vector<Eigen::Vector3d>& pts) const {
    if (pts.size() < 2) return Vector3d(0, 0, 1);
    Vector3d centroid = Vector3d::Zero();
    for (auto& p : pts) centroid += p;
    centroid /= double(pts.size());
    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (auto& p : pts) { Vector3d d = p - centroid; cov += d * d.transpose(); }
    cov /= double(pts.size());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
    return es.eigenvectors().col(2).normalized();
}

// ============================================================================
// 连通分量 (6-邻域)
// ============================================================================
vector<vector<int>> TpmsPipeline::extract_components_6(const VoxelField3D& field) {
    int N = field.total();
    vector<char> visited(N, 0);
    vector<vector<int>> components;
    for (int i = 0; i < N; ++i) {
        if (visited[i] || field.at_idx(i) == 0) continue;
        components.emplace_back();
        queue<int> q;
        q.push(i); visited[i] = 1;
        while (!q.empty()) {
            int cur = q.front(); q.pop();
            components.back().push_back(cur);
            int cx, cy, cz; field.coord(cur, cx, cy, cz);
            for (int d = 0; d < 6; ++d) {
                int nx = cx + k6[d][0], ny = cy + k6[d][1], nz = cz + k6[d][2];
                if (!field.in_bounds(nx, ny, nz)) continue;
                int nid = field.idx(nx, ny, nz);
                if (field.at_idx(nid) != 0 && !visited[nid]) { visited[nid] = 1; q.push(nid); }
            }
        }
    }
    return components;
}

// ============================================================================
// 连通分量 (26-邻域)
// ============================================================================
vector<vector<int>> TpmsPipeline::extract_components_26(
    const vector<uint8_t>& data, int nx, int ny, int nz) const {
    int N = nx * ny * nz;
    vector<char> visited(N, 0);
    vector<vector<int>> components;
    for (int i = 0; i < N; ++i) {
        if (visited[i] || data[i] == 0) continue;
        components.emplace_back();
        queue<int> q; q.push(i); visited[i] = 1;
        while (!q.empty()) {
            int cur = q.front(); q.pop();
            components.back().push_back(cur);
            int cz = cur / (nx * ny), rem = cur % (nx * ny), cy = rem / nx, cx = rem % nx;
            for (int d = 0; d < 26; ++d) {
                int nx2 = cx + k26[d][0], ny2 = cy + k26[d][1], nz2 = cz + k26[d][2];
                if (nx2 < 0 || nx2 >= nx || ny2 < 0 || ny2 >= ny || nz2 < 0 || nz2 >= nz) continue;
                int nid = nx2 + nx * (ny2 + ny * nz2);
                if (data[nid] != 0 && !visited[nid]) { visited[nid] = 1; q.push(nid); }
            }
        }
    }
    return components;
}

// ============================================================================
// 邻居计数
// ============================================================================
int TpmsPipeline::count_neighbors_26(int flat_idx, const VoxelField3D& field) const {
    int cx, cy, cz; field.coord(flat_idx, cx, cy, cz);
    int cnt = 0;
    for (int d = 0; d < 26; ++d) {
        int nx = cx + k26[d][0], ny = cy + k26[d][1], nz = cz + k26[d][2];
        if (field.in_bounds(nx, ny, nz) && field.at(nx, ny, nz) == 1) ++cnt;
    }
    return cnt;
}

int TpmsPipeline::count_neighbors_6(int flat_idx, const VoxelField3D& field) const {
    int cx, cy, cz; field.coord(flat_idx, cx, cy, cz);
    int cnt = 0;
    for (int d = 0; d < 6; ++d) {
        int nx = cx + k6[d][0], ny = cy + k6[d][1], nz = cz + k6[d][2];
        if (field.in_bounds(nx, ny, nz) && field.at(nx, ny, nz) == 1) ++cnt;
    }
    return cnt;
}

// ============================================================
// 3D topology-preserving thinning core: simple point detection
// Uses (26,6) connectivity: foreground 26-connected, background 6-connected
// ============================================================
bool TpmsPipeline::is_simple_point(int flat_idx, const VoxelField3D& field) const {
    int cx, cy, cz;
    field.coord(flat_idx, cx, cy, cz);

    int local_to_flat[27];
    vector<int> fg_locs; fg_locs.reserve(26);

    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                int loc = (dx + 1) + 3 * ((dy + 1) + 3 * (dz + 1));
                if (dx == 0 && dy == 0 && dz == 0) { local_to_flat[loc] = -1; continue; }
                int nx = cx + dx, ny = cy + dy, nz = cz + dz;
                if (field.in_bounds(nx, ny, nz)) {
                    int nid = field.idx(nx, ny, nz);
                    local_to_flat[loc] = nid;
                    if (field.at_idx(nid) == 1) fg_locs.push_back(loc);
                }
                else {
                    local_to_flat[loc] = -2;
                }
            }
        }
    }

    if (fg_locs.empty()) return false;

    // Check 1: foreground 26-connected components == 1
    vector<char> fg_vis(27, 0);
    queue<int> q;
    q.push(fg_locs[0]); fg_vis[fg_locs[0]] = 1;
    int vcnt = 0;
    while (!q.empty()) {
        int cl = q.front(); q.pop(); ++vcnt;
        int clx = cl % 3, cly = (cl / 3) % 3, clz = cl / 9;
        for (int dz = -1; dz <= 1; ++dz)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx) {
                    if (dx == 0 && dy == 0 && dz == 0) continue;
                    int nlx = clx + dx, nly = cly + dy, nlz = clz + dz;
                    if (nlx < 0 || nlx > 2 || nly < 0 || nly > 2 || nlz < 0 || nlz > 2) continue;
                    int nl = nlx + 3 * (nly + 3 * nlz);
                    if (fg_vis[nl]) continue;
                    int nid = local_to_flat[nl];
                    if (nid >= 0 && field.at_idx(nid) == 1) { fg_vis[nl] = 1; q.push(nl); }
                }
    }
    if (vcnt != (int)fg_locs.size()) return false;

    // 检查2: 与 p 6-相邻的背景点 在 3x3x3 中 6-连通分量 = 1
    vector<int> bg_6adj;
    for (int d = 0; d < 6; ++d) {
        int nl = (1 + k6[d][0]) + 3 * ((1 + k6[d][1]) + 3 * (1 + k6[d][2]));
        int nid = local_to_flat[nl];
        if (nid == -2 || (nid >= 0 && field.at_idx(nid) == 0)) bg_6adj.push_back(nl);
    }
    if (bg_6adj.empty()) return false;

    vector<char> bg_vis(27, 0);
    q.push(bg_6adj[0]); bg_vis[bg_6adj[0]] = 1;
    while (!q.empty()) {
        int cl = q.front(); q.pop();
        int clx = cl % 3, cly = (cl / 3) % 3, clz = cl / 9;
        for (int d = 0; d < 6; ++d) {
            int nlx = clx + k6[d][0], nly = cly + k6[d][1], nlz = clz + k6[d][2];
            if (nlx < 0 || nlx > 2 || nly < 0 || nly > 2 || nlz < 0 || nlz > 2) continue;
            int nl = nlx + 3 * (nly + 3 * nlz);
            if (bg_vis[nl] || nl == 13) continue;
            int nid = local_to_flat[nl];
            bool bg = (nid == -2) || (nid >= 0 && field.at_idx(nid) == 0);
            if (bg) { bg_vis[nl] = 1; q.push(nl); }
        }
    }
    for (int loc : bg_6adj) if (!bg_vis[loc]) return false;
    return true;
}

bool TpmsPipeline::is_endpoint(int flat_idx) const {
    int N = params_.resolution;
    int cz = flat_idx / (N * N), rem = flat_idx % (N * N), cy = rem / N, cx = rem % N;
    int cnt = 0;
    for (int d = 0; d < 26; ++d) {
        int nx = cx + k26[d][0], ny = cy + k26[d][1], nz = cz + k26[d][2];
        if (nx < 0 || nx >= N || ny < 0 || ny >= N || nz < 0 || nz >= N) continue;
        if (skeleton_data_[nx + N * (ny + N * nz)] == 1) ++cnt;
    }
    return cnt == 1;
}

void TpmsPipeline::distance_transform(const VoxelField3D& src, vector<double>& dist) const {
    int N = src.total(); dist.assign(N, 0.0);
    queue<int> q; vector<char> visited(N, 0);
    for (int i = 0; i < N; ++i) {
        if (src.at_idx(i) == 0) continue;
        int cx, cy, cz; src.coord(i, cx, cy, cz);
        bool bdry = false;
        for (int d = 0; d < 6; ++d) {
            int nx = cx + k6[d][0], ny = cy + k6[d][1], nz = cz + k6[d][2];
            if (!src.in_bounds(nx, ny, nz) || src.at(nx, ny, nz) == 0) { bdry = true; break; }
        }
        if (bdry) { dist[i] = 1.0; q.push(i); visited[i] = 1; }
    }
    while (!q.empty()) {
        int cur = q.front(); q.pop();
        int cx, cy, cz; src.coord(cur, cx, cy, cz);
        for (int d = 0; d < 6; ++d) {
            int nx = cx + k6[d][0], ny = cy + k6[d][1], nz = cz + k6[d][2];
            if (!src.in_bounds(nx, ny, nz)) continue;
            int nid = src.idx(nx, ny, nz);
            if (src.at_idx(nid) == 0 || visited[nid]) continue;
            dist[nid] = dist[cur] + 1.0; visited[nid] = 1; q.push(nid);
        }
    }
}

// ============================================================================
// Step 1: 生成 TPMS 场 (Gyroid)
// ============================================================================
void TpmsPipeline::step1_generateTpmsField() {
    const int R = params_.resolution;
    const double sp = params_.spacing();
    const double ox = -params_.half_bbox;
    const int total = params_.grid_total();
    tpms_field_.resize(total);

    double ax = params_.alpha_x, ay = params_.alpha_y, az = params_.alpha_z;
    double px = params_.phase_x, py = params_.phase_y, pz = params_.phase_z;

    cout << "[Step 1] Generating ";
    if (params_.is_perturbed) cout << "PERTURBED ";
    cout << "Gyroid TPMS at " << R << "^3" << endl;
    if (params_.is_perturbed)
        cout << "  alpha=(" << ax << "," << ay << "," << az
        << ") phase=(" << px << "," << py << "," << pz << ")" << endl;

    for (int z = 0; z < R; ++z) {
        double wz = az * (ox + z * sp) + pz;
        for (int y = 0; y < R; ++y) {
            double wy = ay * (ox + y * sp) + py;
            for (int x = 0; x < R; ++x) {
                double wx = ax * (ox + x * sp) + px;
                double f = sin(wx) * cos(wy) + sin(wy) * cos(wz) + sin(wz) * cos(wx);
                tpms_field_[x + R * (y + R * z)] = f;
            }
        }
    }

    auto [mn, mx] = minmax_element(tpms_field_.begin(), tpms_field_.end());
    cout << "  Field range: [" << *mn << ", " << *mx << "]" << endl;
}

// ============================================================================
// Step 2: 体素化
// ============================================================================
void TpmsPipeline::step2_voxelize() {
    const int R = params_.resolution;
    const double margin = params_.spacing() * 0.5;
    void_mask_.nx = R; void_mask_.ny = R; void_mask_.nz = R;
    void_mask_.data.resize(params_.grid_total(), 0);
    int vc = 0;
    for (int i = 0; i < params_.grid_total(); ++i) {
        if (tpms_field_[i] < params_.iso_value - margin) { void_mask_.data[i] = 1; ++vc; }
    }
    stats_.total_void_voxels = vc;
    cout << "[Step 2] Voxelized: " << vc << " void / " << params_.grid_total()
        << " (" << 100.0 * vc / params_.grid_total() << "%)" << endl;
}

// ============================================================================
// Step 3: 提取最大连通 void region
// ============================================================================
void TpmsPipeline::step3_extractLargestVoidRegion() {
    cout << "[Step 3] Extracting largest connected void region..." << endl;
    auto comps = extract_components_6(void_mask_);
    if (comps.empty()) { cerr << "  ERROR: No void components!" << endl; return; }
    cout << "  Found " << comps.size() << " components." << endl;
    size_t mi = 0;
    for (size_t i = 1; i < comps.size(); ++i)
        if (comps[i].size() > comps[mi].size()) mi = i;
    largest_void_.nx = params_.resolution; largest_void_.ny = params_.resolution; largest_void_.nz = params_.resolution;
    largest_void_.data.resize(params_.grid_total(), 0);
    for (int idx : comps[mi]) largest_void_.data[idx] = 1;
    stats_.largest_void_voxels = (int)comps[mi].size();
    cout << "  Largest: " << stats_.largest_void_voxels << " voxels" << endl;
    vector<size_t> cs; for (auto& c : comps) cs.push_back(c.size());
    sort(cs.begin(), cs.end(), greater<size_t>());
    cout << "  Top sizes: "; for (size_t i = 0; i < min(cs.size(), size_t(5)); ++i) cout << cs[i] << " ";
    cout << endl;
}

// ============================================================================
// Step 4: 3D 骨架提取 (Topology-Preserving Thinning)
// ============================================================================
void TpmsPipeline::step4_skeletonize() {
    const int R = params_.resolution, total = params_.grid_total();
    cout << "[Step 4] 3D skeletonization..." << endl;

    skeleton_data_.resize(total);
    for (int i = 0; i < total; ++i) skeleton_data_[i] = largest_void_.data[i];

    auto mkview = [&]() { VoxelField3D v; v.nx = R; v.ny = R; v.nz = R; v.data = skeleton_data_; return v; };

    static const int dirs[6][3] = { {0,-1,0},{0,1,0},{1,0,0},{-1,0,0},{0,0,1},{0,0,-1} };
    int total_removed = 0;

    for (int it = 0; it < R; ++it) {
        int rem_iter = 0;
        for (int di = 0; di < 6; ++di) {
            int dx = dirs[di][0], dy = dirs[di][1], dz = dirs[di][2];
            vector<int> todel;
            for (int i = 0; i < total; ++i) {
                if (skeleton_data_[i] == 0) continue;
                int cx, cy, cz; void_mask_.coord(i, cx, cy, cz);
                int nx = cx + dx, ny = cy + dy, nz = cz + dz;
                if (nx >= 0 && nx < R && ny >= 0 && ny < R && nz >= 0 && nz < R) {
                    if (skeleton_data_[nx + R * (ny + R * nz)] == 1) continue;
                }
                if (is_endpoint(i)) continue;
                VoxelField3D vw = mkview();
                if (is_simple_point(i, vw)) todel.push_back(i);
            }
            for (int idx : todel) { skeleton_data_[idx] = 0; ++rem_iter; }
        }
        total_removed += rem_iter;
        if (rem_iter == 0) { cout << "  Converged at iter " << it + 1 << endl; break; }
        if ((it + 1) % 10 == 0)
            cout << "  Iter " << it + 1 << ": -" << rem_iter << " (total -" << total_removed << ")" << endl;
    }

    int sk = 0;
    for (int i = 0; i < total; ++i) if (skeleton_data_[i] == 1) ++sk;
    stats_.skeleton_voxels_initial = sk;

    double sp = params_.spacing(), ox = -params_.half_bbox;
    skeleton_voxels_.clear(); skeleton_voxels_.reserve(sk);
    for (int i = 0; i < total; ++i) {
        if (skeleton_data_[i] == 0) continue;
        SkeletonVoxel sv; sv.idx_flat = i;
        void_mask_.coord(i, sv.x, sv.y, sv.z);
        sv.world_pos = Vector3d(ox + sv.x * sp, ox + sv.y * sp, ox + sv.z * sp);
        int deg = 0;
        for (int d = 0; d < 26; ++d) {
            int nx = sv.x + k26[d][0], ny = sv.y + k26[d][1], nz = sv.z + k26[d][2];
            if (nx < 0 || nx >= R || ny < 0 || ny >= R || nz < 0 || nz >= R) continue;
            if (skeleton_data_[nx + R * (ny + R * nz)] == 1) ++deg;
        }
        sv.degree = deg; sv.dist_to_surface = 0.0;
        skeleton_voxels_.push_back(sv);
    }
    cout << "  Done: " << sk << " skeleton voxels (removed " << total_removed << ")" << endl;
}

// ============================================================================
// Step 5: 清理短枝
// ============================================================================
void TpmsPipeline::step5_pruneSkeleton(double min_branch_len) {
    const int R = params_.resolution, total = params_.grid_total();
    cout << "[Step 5] Pruning short branches (min=" << min_branch_len << ")..." << endl;
    int pruned = 0;
    for (int it = 0; it < 30; ++it) {
        vector<int> eps;
        for (int i = 0; i < total; ++i) if (skeleton_data_[i] == 1 && is_endpoint(i)) eps.push_back(i);
        if (eps.empty()) break;
        int it_pruned = 0;
        for (int ep : eps) {
            if (skeleton_data_[ep] == 0) continue;
            vector<int> branch; branch.push_back(ep);
            int cur = ep; set<int> vis; vis.insert(cur);
            bool to_junc = false;
            while (true) {
                int cx, cy, cz; void_mask_.coord(cur, cx, cy, cz);
                bool moved = false;
                for (int d = 0; d < 26; ++d) {
                    int nx = cx + k26[d][0], ny = cy + k26[d][1], nz = cz + k26[d][2];
                    if (nx < 0 || nx >= R || ny < 0 || ny >= R || nz < 0 || nz >= R) continue;
                    int nid = nx + R * (ny + R * nz);
                    if (skeleton_data_[nid] == 0 || vis.count(nid)) continue;
                    cur = nid; moved = true; break;
                }
                if (!moved) break;
                vis.insert(cur); branch.push_back(cur);
                int deg = 0;
                int dcx, dcy, dcz; void_mask_.coord(cur, dcx, dcy, dcz);
                for (int d = 0; d < 26; ++d) {
                    int nx = dcx + k26[d][0], ny = dcy + k26[d][1], nz = dcz + k26[d][2];
                    if (nx < 0 || nx >= R || ny < 0 || ny >= R || nz < 0 || nz >= R) continue;
                    if (skeleton_data_[nx + R * (ny + R * nz)] == 1) ++deg;
                }
                if (deg >= 3) { to_junc = true; break; }
                if (deg == 1 && branch.size() > 1) break;
                if ((int)branch.size() > R) break;
            }
            int last_deg = 0;
            if (!branch.empty()) {
                int lx, ly, lz; void_mask_.coord(branch.back(), lx, ly, lz);
                for (int d = 0; d < 26; ++d) {
                    int nx = lx + k26[d][0], ny = ly + k26[d][1], nz = lz + k26[d][2];
                    if (nx < 0 || nx >= R || ny < 0 || ny >= R || nz < 0 || nz >= R) continue;
                    if (skeleton_data_[nx + R * (ny + R * nz)] == 1) ++last_deg;
                }
            }
            if (to_junc && (int)branch.size() <= (int)min_branch_len) {
                for (int i = 0; i < (int)branch.size() - 1; ++i) {
                    if (skeleton_data_[branch[i]] == 1) { skeleton_data_[branch[i]] = 0; ++it_pruned; }
                }
            }
        }
        pruned += it_pruned;
        if (it_pruned == 0) break;
    }
    cout << "  Pruned " << pruned << " voxels total" << endl;

    int sk = 0; double sp = params_.spacing(), ox = -params_.half_bbox;
    for (int i = 0; i < total; ++i) if (skeleton_data_[i] == 1) ++sk;
    stats_.skeleton_voxels_pruned = sk;
    skeleton_voxels_.clear(); skeleton_voxels_.reserve(sk);
    for (int i = 0; i < total; ++i) {
        if (skeleton_data_[i] == 0) continue;
        SkeletonVoxel sv; sv.idx_flat = i;
        void_mask_.coord(i, sv.x, sv.y, sv.z);
        sv.world_pos = Vector3d(ox + sv.x * sp, ox + sv.y * sp, ox + sv.z * sp);
        int deg = 0;
        for (int d = 0; d < 26; ++d) {
            int nx = sv.x + k26[d][0], ny = sv.y + k26[d][1], nz = sv.z + k26[d][2];
            if (nx < 0 || nx >= R || ny < 0 || ny >= R || nz < 0 || nz >= R) continue;
            if (skeleton_data_[nx + R * (ny + R * nz)] == 1) ++deg;
        }
        sv.degree = deg; sv.dist_to_surface = 0.0;
        skeleton_voxels_.push_back(sv);
    }
    cout << "  " << sk << " skeleton voxels after pruning" << endl;
}

// ============================================================================
// Step 6: Skeleton → Graph
// ============================================================================
void TpmsPipeline::step6_buildGraph(double junc_radius, double max_edge_no_node) {
    const int R = params_.resolution, total = params_.grid_total();
    const double sp = params_.spacing(), ox = -params_.half_bbox;
    cout << "[Step 6] Building graph from skeleton..." << endl;

    // 6.1 分类每个骨架体素
    vector<int> vtype(total, -1);
    vector<int> eps_flat, juncs_flat;
    for (auto& sv : skeleton_voxels_) {
        int i = sv.idx_flat;
        if (sv.degree == 1) { vtype[i] = 1; eps_flat.push_back(i); }
        else if (sv.degree >= 3) { vtype[i] = 2; juncs_flat.push_back(i); }
        else vtype[i] = 0;
    }
    cout << "  Endpoints: " << eps_flat.size() << ", Junction voxels: " << juncs_flat.size() << endl;

    // 6.2 Junction clustering (26-连通)
    vector<vector<int>> jclusters;
    {
        vector<char> jv(total, 0);
        for (int j : juncs_flat) {
            if (jv[j]) continue;
            jclusters.emplace_back();
            queue<int> q; q.push(j); jv[j] = 1;
            while (!q.empty()) {
                int cur = q.front(); q.pop();
                jclusters.back().push_back(cur);
                int cx, cy, cz; void_mask_.coord(cur, cx, cy, cz);
                for (int d = 0; d < 26; ++d) {
                    int nx = cx + k26[d][0], ny = cy + k26[d][1], nz = cz + k26[d][2];
                    if (nx < 0 || nx >= R || ny < 0 || ny >= R || nz < 0 || nz >= R) continue;
                    int nid = nx + R * (ny + R * nz);
                    if (vtype[nid] >= 2 && !jv[nid]) { jv[nid] = 1; q.push(nid); }
                }
            }
        }
    }
    cout << "  Junction clusters: " << jclusters.size() << endl;

    // 6.3 创建 graph nodes
    vector<int> node_of_flat(total, -1);
    graph_nodes_.clear();
    graph_nodes_.reserve(jclusters.size() + eps_flat.size() + 50);

    for (auto& cl : jclusters) {
        GraphNode gn; gn.id = (int)graph_nodes_.size(); gn.type = GraphNodeType::JUNCTION;
        Vector3d c = Vector3d::Zero();
        for (int f : cl) { int cx, cy, cz; void_mask_.coord(f, cx, cy, cz); c += Vector3d(ox + cx * sp, ox + cy * sp, ox + cz * sp); }
        gn.position = c / double(cl.size());
        for (int f : cl) node_of_flat[f] = gn.id;
        graph_nodes_.push_back(gn);
    }
    for (int ep : eps_flat) {
        if (node_of_flat[ep] >= 0) continue;
        GraphNode gn; gn.id = (int)graph_nodes_.size(); gn.type = GraphNodeType::ENDPOINT;
        int cx, cy, cz; void_mask_.coord(ep, cx, cy, cz);
        gn.position = Vector3d(ox + cx * sp, ox + cy * sp, ox + cz * sp);
        node_of_flat[ep] = gn.id;
        graph_nodes_.push_back(gn);
    }
    cout << "  Total graph nodes: " << graph_nodes_.size() << endl;

    // 6.4 追踪骨架边
    graph_edges_.clear();
    set<pair<int, int>> edge_dedup;

    for (int ni = 0; ni < (int)graph_nodes_.size(); ++ni) {
        auto& node = graph_nodes_[ni];
        int seed = -1;
        for (auto& cl : jclusters) {
            for (int f : cl) if (node_of_flat[f] == node.id) { seed = f; break; }
            if (seed >= 0) break;
        }
        if (seed < 0) {
            for (int ep : eps_flat) if (node_of_flat[ep] == node.id) { seed = ep; break; }
        }
        if (seed < 0) continue;

        int sx, sy, sz; void_mask_.coord(seed, sx, sy, sz);
        for (int d = 0; d < 26; ++d) {
            int nx = sx + k26[d][0], ny = sy + k26[d][1], nz = sz + k26[d][2];
            if (nx < 0 || nx >= R || ny < 0 || ny >= R || nz < 0 || nz >= R) continue;
            int nid = nx + R * (ny + R * nz);
            if (skeleton_data_[nid] == 0) continue;
            if (node_of_flat[nid] >= 0) {
                int an = node_of_flat[nid];
                pair<int, int> key(min(ni, an), max(ni, an));
                if (ni != an && !edge_dedup.count(key)) {
                    edge_dedup.insert(key);
                    GraphEdge ge; ge.id = (int)graph_edges_.size();
                    ge.node_a = ni; ge.node_b = an;
                    ge.polyline.push_back(graph_nodes_[ni].position);
                    ge.polyline.push_back(graph_nodes_[an].position);
                    ge.length_world = (ge.polyline[1] - ge.polyline[0]).norm();
                    graph_edges_.push_back(ge);
                }
                continue;
            }
            // 追踪路径 path
            vector<int> path; path.push_back(nid);
            int cur = nid; set<int> pvis; pvis.insert(cur); pvis.insert(seed);
            bool reached = false; int target_node = -1;
            while (true) {
                int cx, cy, cz; void_mask_.coord(cur, cx, cy, cz);
                int next = -1;
                for (int dd = 0; dd < 26; ++dd) {
                    int nnx = cx + k26[dd][0], nny = cy + k26[dd][1], nnz = cz + k26[dd][2];
                    if (nnx < 0 || nnx >= R || nny < 0 || nny >= R || nnz < 0 || nnz >= R) continue;
                    int nnid = nnx + R * (nny + R * nnz);
                    if (skeleton_data_[nnid] == 0 || pvis.count(nnid)) continue;
                    if (node_of_flat[nnid] >= 0) { target_node = node_of_flat[nnid]; reached = true; break; }
                    next = nnid; break;
                }
                if (reached) break;
                if (next < 0) break;
                cur = next; pvis.insert(cur); path.push_back(cur);
                if ((int)path.size() > R * 3) break;
            }
            if (!reached || target_node < 0) continue;
            pair<int, int> key(min(ni, target_node), max(ni, target_node));
            if (ni == target_node || edge_dedup.count(key)) continue;
            edge_dedup.insert(key);
            GraphEdge ge; ge.id = (int)graph_edges_.size();
            ge.node_a = ni; ge.node_b = target_node;
            ge.polyline.push_back(graph_nodes_[ni].position);
            for (int f : path) {
                int px, py, pz; void_mask_.coord(f, px, py, pz);
                ge.polyline.push_back(Vector3d(ox + px * sp, ox + py * sp, ox + pz * sp));
            }
            ge.polyline.push_back(graph_nodes_[target_node].position);
            ge.length_world = 0;
            for (size_t i = 1; i < ge.polyline.size(); ++i)
                ge.length_world += (ge.polyline[i] - ge.polyline[i - 1]).norm();
            graph_edges_.push_back(ge);
        }
    }

    // 6.5 长边插入 virtual nodes
    int vn_added = 0;
    vector<GraphEdge> new_edges;
    for (auto& e : graph_edges_) {
        if (e.length_world > max_edge_no_node * sp && e.polyline.size() > 2) {
            double seg = e.length_world / ceil(e.length_world / (max_edge_no_node * sp));
            vector<int> inserted;
            double acc = 0; size_t si = 1;
            while (acc + seg < e.length_world && si < e.polyline.size()) {
                double seg_len = (e.polyline[si] - e.polyline[si - 1]).norm();
                while (acc + seg_len >= seg && acc + seg < e.length_world && si < e.polyline.size()) {
                    double t = (seg - acc) / seg_len;
                    Vector3d vp = e.polyline[si - 1] + t * (e.polyline[si] - e.polyline[si - 1]);
                    GraphNode vn; vn.id = (int)graph_nodes_.size(); vn.type = GraphNodeType::VIRTUAL;
                    vn.position = vp;
                    graph_nodes_.push_back(vn);
                    inserted.push_back(vn.id); ++vn_added;
                    acc += seg; seg = e.length_world / ceil(e.length_world / (max_edge_no_node * sp));
                }
                acc += seg_len; ++si;
            }
            int prev = e.node_a;
            for (int vid : inserted) {
                GraphEdge ne; ne.id = (int)graph_edges_.size();
                ne.node_a = prev; ne.node_b = vid;
                ne.polyline.push_back(graph_nodes_[prev].position);
                ne.polyline.push_back(graph_nodes_[vid].position);
                ne.length_world = (ne.polyline[1] - ne.polyline[0]).norm();
                ne.id = (int)new_edges.size(); new_edges.push_back(ne);
                prev = vid;
            }
            GraphEdge ne; ne.id = (int)graph_edges_.size();
            ne.node_a = prev; ne.node_b = e.node_b;
            ne.polyline.push_back(graph_nodes_[prev].position);
            ne.polyline.push_back(graph_nodes_[e.node_b].position);
            ne.length_world = (ne.polyline[1] - ne.polyline[0]).norm();
            new_edges.push_back(ne);
        }
        else {
            new_edges.push_back(e);
        }
    }
    graph_edges_ = move(new_edges);

    // 更新 node neighbors
    for (auto& n : graph_nodes_) n.neighbors.clear();
    for (int ei = 0; ei < (int)graph_edges_.size(); ++ei) {
        auto& e = graph_edges_[ei]; e.id = ei;
        graph_nodes_[e.node_a].neighbors.push_back(ei);
        graph_nodes_[e.node_b].neighbors.push_back(ei);
    }
    cout << "  Virtual nodes added: " << vn_added << endl;
    cout << "  Final graph: " << graph_nodes_.size() << " nodes, " << graph_edges_.size() << " edges" << endl;
}

// ============================================================================
// Step 7: VP 计算
// ============================================================================
void TpmsPipeline::step7_computeVP() {
    const int N = (int)graph_nodes_.size();
    const int R = params_.resolution;
    const double sp = params_.spacing(), ox = -params_.half_bbox;
    cout << "[Step 7] Computing Visual Permeability (VP) on graph..." << endl;

    // 7.0 Node location classification (surface vs interior)
    // distance from each node to nearest TPMS wall
    const double surf_threshold = 2.0 * sp; // ~2 voxels

    for (auto& node : graph_nodes_) {
        // sample SDF at node position
        double min_abs_f = 1e9;
        Vector3d p = node.position;
        // interpolate on voxel grid
        double fx = (p.x() - ox) / sp;
        double fy = (p.y() - ox) / sp;
        double fz = (p.z() - ox) / sp;
        int ix = (int)round(fx), iy = (int)round(fy), iz = (int)round(fz);
        ix = clamp(ix, 0, R - 1); iy = clamp(iy, 0, R - 1); iz = clamp(iz, 0, R - 1);
        int fidx = ix + R * (iy + R * iz);
        double fval = abs(tpms_field_[fidx]);
        node.dist_to_surface = fval; // |f| approximates distance to isosurface
        if (fval < surf_threshold)
            node.location = NodeLocation::SURFACE;
        else
            node.location = NodeLocation::INTERIOR;
    }

    int surfcnt = 0;
    for (auto& n : graph_nodes_) if (n.location == NodeLocation::SURFACE) ++surfcnt;
    stats_.surface_node_count = surfcnt;
    stats_.interior_node_count = N - surfcnt;
    cout << "  Surface nodes: " << surfcnt << ", Interior nodes: " << (N - surfcnt) << endl;

    // 7.1 构建邻接表
    vector<vector<int>> adj(N);
    for (auto& e : graph_edges_) {
        adj[e.node_a].push_back(e.node_b);
        adj[e.node_b].push_back(e.node_a);
    }

    vector<bool> is_surface(N, false);
    for (auto& n : graph_nodes_) if (n.location == NodeLocation::SURFACE) is_surface[n.id] = true;

    // 7.2 对每个节点计算 VP
    vp_result_.node_vp.assign(N, 0.0);
    vp_result_.node_best_path.resize(N);
    double total_vp = 0.0;

    // 策略: 对每个节点, BFS 找到可达的 surface nodes, 选择 pair 使路径 VP 最大
    for (int ni = 0; ni < N; ++ni) {
        // BFS from ni to collect surface nodes and distances
        vector<int> dist(N, -1), parent(N, -1);
        queue<int> q;
        dist[ni] = 0; q.push(ni);

        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : adj[u]) {
                if (dist[v] >= 0) continue;
                dist[v] = dist[u] + 1; parent[v] = u; q.push(v);
            }
        }

        // 收集可达的 surface nodes
        vector<int> surf_nodes;
        for (int j = 0; j < N; ++j) {
            if (is_surface[j] && dist[j] > 0) surf_nodes.push_back(j);
        }
        if (surf_nodes.size() < 2) {
            vp_result_.node_vp[ni] = 0.0;
            continue;
        }

        // 启发式: 选距离最远的两个 surface node
        int s1 = surf_nodes[0], s2 = surf_nodes[0];
        for (int s : surf_nodes) {
            if (dist[s] > dist[s1]) s1 = s;
        }
        // BFS from s1 找最远的 surface node
        vector<int> dist2(N, -1), parent2(N, -1);
        queue<int> q2; dist2[s1] = 0; q2.push(s1);
        while (!q2.empty()) { int u = q2.front(); q2.pop(); for (int v : adj[u]) if (dist2[v] < 0) { dist2[v] = dist2[u] + 1; parent2[v] = u; q2.push(v); } }
        for (int s : surf_nodes) { if (dist2[s] > dist2[s2]) s2 = s; }

        // 重建 s1→ni→s2 路径
        vector<int> path;
        for (int cur = ni; cur != s1; cur = parent[cur]) path.push_back(cur);
        path.push_back(s1);
        reverse(path.begin(), path.end());
        for (int cur = parent2[ni]; cur != s2; cur = parent2[cur]) path.push_back(cur);
        path.push_back(s2);

        // 去重
        vector<int> dedup;
        for (int v : path) { if (dedup.empty() || dedup.back() != v) dedup.push_back(v); }

        // 计算 VP
        vector<Vector3d> path_pts;
        for (int v : dedup) path_pts.push_back(graph_nodes_[v].position);
        double score = compute_path_vp(path_pts, is_surface);
        vp_result_.node_vp[ni] = score;
        vp_result_.node_best_path[ni] = dedup;
        total_vp += score;
    }

    vp_result_.total_vp = total_vp / double(N);
    stats_.final_vp = vp_result_.total_vp;

    // compute average degree and stats
    double sum_deg = 0;
    for (auto& node : graph_nodes_) sum_deg += (double)node.neighbors.size();
    stats_.avg_degree = sum_deg / double(N);
    stats_.graph_node_count = N;
    stats_.graph_edge_count = (int)graph_edges_.size();

    std::cout << "  Total VP: " << stats_.final_vp << std::endl;
    std::cout << "  Avg degree: " << stats_.avg_degree << std::endl;
}

double TpmsPipeline::compute_path_vp(
    const vector<Vector3d>& path_pts,
    const vector<bool>& is_surface) const
{
    int psize = (int)path_pts.size();
    if (psize < 2) return 0.0;

    // T_angle
    double T_angle;
    const double angle_floor = 0.6;
    if (psize == 2) {
        T_angle = angle_floor;
    }
    else {
        double prod = 1.0;
        for (int i = 1; i < psize - 1; ++i) {
            double ad = abs_angle_deg(path_pts[i - 1] - path_pts[i], path_pts[i + 1] - path_pts[i]) / 180.0;
            prod *= ad;
        }
        T_angle = pow(prod, 0.5);
    }

    // T_length
    double L0 = 5.0;
    double T_length = 1.0 - exp(-double(psize - 1) / L0);

    // T_location
    int inner_cnt = 0;
    for (int i = 1; i < psize - 1; ++i) {
        // 这里简化: 平均节点位置的 SDF
    }
    // 使用路径点与 tpms_field_ 的近似
    int surf_line_cnt = 0;
    const double loc_thres = 2.0 * params_.spacing();
    for (size_t i = 0; i < path_pts.size() - 1; ++i) {
        Vector3d mid = 0.5 * (path_pts[i] + path_pts[i + 1]);
        double fx = (mid.x() + params_.half_bbox) / params_.spacing();
        double fy = (mid.y() + params_.half_bbox) / params_.spacing();
        double fz = (mid.z() + params_.half_bbox) / params_.spacing();
        int ix = clamp((int)fx, 0, params_.resolution - 1);
        int iy = clamp((int)fy, 0, params_.resolution - 1);
        int iz = clamp((int)fz, 0, params_.resolution - 1);
        double fv = abs(tpms_field_[ix + params_.resolution * (iy + params_.resolution * iz)]);
        if (fv >= loc_thres) ++surf_line_cnt;
    }
    int inner_line_cnt = max(0, psize - 1 - surf_line_cnt);
    double r_inner = (psize == 2) ? 0 : double(inner_line_cnt) / double(psize - 1);
    double r_surface = (psize == 2) ? 0 : double(surf_line_cnt) / double(psize - 1);
    double T_location = 1.0 / (1.0 + exp(-8.0 * (r_inner - r_surface)));

    // S_horiz
    Vector3d dir = computePrincipalDirection(path_pts);
    double S_horiz = 1.0 - abs(dir.dot(Vector3d(0, 0, 1)));

    double score = vp_weights_.w_angle * T_angle
        + vp_weights_.w_length * T_length
        + vp_weights_.w_location * T_location
        + vp_weights_.w_direction * S_horiz;
    return score;
}

// ============================================================================
// Step 8: export visualization data
// ============================================================================
void TpmsPipeline::step8_exportData(const string& output_prefix) {
    std::cout << "[Step 8] Exporting visualization data to: " << output_prefix << std::endl;

    {   // 8.1 TPMS field -> NPY
        string fn = output_prefix + "_tpms_field.npy";
        ofstream f(fn, ios::binary);
        if (!f) { std::cerr << "Failed to open " << fn << std::endl; return; }
        f.write("\x93NUMPY", 6);
        f.write("\x01\x00", 2);
        int R = params_.resolution;
        string shape = "(" + to_string(R) + "," + to_string(R) + "," + to_string(R) + ")";
        string header = "{'descr':'<f8','fortran_order':False,'shape':" + shape + "}";
        int pad = (16 - ((10 + (int)header.size()) % 16)) % 16;
        header.append(pad, ' ');
        header += "\n";
        f.write(header.data(), header.size());
        f.write((const char*)tpms_field_.data(), tpms_field_.size() * sizeof(double));
        std::cout << "  Exported: " << fn << " (" << tpms_field_.size() << " doubles)" << std::endl;
    }

    {   // 8.2 Skeleton voxels -> CSV
        string fn = output_prefix + "_skeleton.csv";
        ofstream f(fn);
        if (!f) { std::cerr << "Failed to open " << fn << std::endl; return; }
        f << "x,y,z,degree,dist_to_surface\n";
        for (auto& sv : skeleton_voxels_) {
            f << sv.world_pos.x() << "," << sv.world_pos.y() << "," << sv.world_pos.z()
                << "," << sv.degree << "," << sv.dist_to_surface << "\n";
        }
        std::cout << "  Exported: " << fn << " (" << skeleton_voxels_.size() << " points)" << std::endl;
    }

    {   // 8.3 Skeleton mask -> NPY
        string fn = output_prefix + "_skeleton.npy";
        ofstream f(fn, ios::binary);
        if (!f) { std::cerr << "Failed to open " << fn << std::endl; return; }
        f.write("\x93NUMPY", 6);
        f.write("\x01\x00", 2);
        int R = params_.resolution;
        string shape = "(" + to_string(R) + "," + to_string(R) + "," + to_string(R) + ")";
        string header = "{'descr':'<u1','fortran_order':False,'shape':" + shape + "}";
        int pad = (16 - ((10 + (int)header.size()) % 16)) % 16;
        header.append(pad, ' ');
        header += "\n";
        f.write(header.data(), header.size());
        f.write((const char*)skeleton_data_.data(), skeleton_data_.size());
        std::cout << "  Exported: " << fn << " (" << skeleton_data_.size() << " bytes)" << std::endl;
    }

    {   // 8.4 Graph nodes -> OBJ
        string fn = output_prefix + "_graph_nodes.obj";
        ofstream f(fn);
        if (!f) { std::cerr << "Failed to open " << fn << std::endl; return; }
        f << "# Graph nodes: " << graph_nodes_.size() << "\n";
        f << "# surface=" << stats_.surface_node_count << " interior=" << stats_.interior_node_count << "\n";
        for (auto& nd : graph_nodes_) {
            f << "v " << nd.position.x() << " " << nd.position.y() << " " << nd.position.z()
                << " # id=" << nd.id << " type=" << (int)nd.type
                << " loc=" << (nd.location == NodeLocation::SURFACE ? "SURFACE" : "INTERIOR")
                << " dist=" << nd.dist_to_surface << "\n";
        }
        std::cout << "  Exported: " << fn << " (" << graph_nodes_.size() << " vertices)" << std::endl;
    }

    {   // 8.5 Graph edges -> OBJ (lines)
        string fn = output_prefix + "_graph_edges.obj";
        ofstream f(fn);
        if (!f) { std::cerr << "Failed to open " << fn << std::endl; return; }
        f << "# Graph edges: " << graph_edges_.size() << "\n";
        for (auto& nd : graph_nodes_) {
            f << "v " << nd.position.x() << " " << nd.position.y() << " " << nd.position.z() << "\n";
        }
        for (auto& e : graph_edges_) {
            f << "l " << (e.node_a + 1) << " " << (e.node_b + 1) << "\n";
        }
        std::cout << "  Exported: " << fn << " (" << graph_edges_.size() << " lines)" << std::endl;
    }

    {   // 8.6 Graph edges with polyline waypoints -> OBJ
        string fn = output_prefix + "_graph_polylines.obj";
        ofstream f(fn);
        if (!f) { std::cerr << "Failed to open " << fn << std::endl; return; }
        f << "# Graph polylines\n";
        int voff = 0;
        for (auto& e : graph_edges_) {
            for (auto& pt : e.polyline) {
                f << "v " << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
            }
            for (size_t i = 0; i + 1 < e.polyline.size(); ++i) {
                f << "l " << (voff + (int)i + 1) << " " << (voff + (int)i + 2) << "\n";
            }
            voff += (int)e.polyline.size();
        }
        std::cout << "  Exported: " << fn << " (" << voff << " polyline vertices)" << std::endl;
    }

    {   // 8.7 Largest void region -> NPY
        string fn = output_prefix + "_void_region.npy";
        ofstream f(fn, ios::binary);
        if (!f) { std::cerr << "Failed to open " << fn << std::endl; return; }
        f.write("\x93NUMPY", 6);
        f.write("\x01\x00", 2);
        int R = params_.resolution;
        string shape = "(" + to_string(R) + "," + to_string(R) + "," + to_string(R) + ")";
        string header = "{'descr':'<u1','fortran_order':False,'shape':" + shape + "}";
        int pad = (16 - ((10 + (int)header.size()) % 16)) % 16;
        header.append(pad, ' ');
        header += "\n";
        f.write(header.data(), header.size());
        f.write((const char*)largest_void_.data.data(), largest_void_.data.size());
        std::cout << "  Exported: " << fn << " (" << largest_void_.data.size() << " bytes)" << std::endl;
    }

    {   // 8.8 Stats summary
        string fn = output_prefix + "_stats.txt";
        ofstream f(fn);
        if (!f) { std::cerr << "Failed to open " << fn << std::endl; return; }
        f << "# TPMS Pipeline Statistics\n";
        f << "perturbed=" << (stats_.is_perturbed ? "yes" : "no") << "\n";
        f << "resolution=" << stats_.resolution << "\n";
        f << "bbox_half=" << stats_.bbox_half << "\n";
        f << "total_void_voxels=" << stats_.total_void_voxels << "\n";
        f << "largest_void_voxels=" << stats_.largest_void_voxels << "\n";
        f << "skeleton_voxels_initial=" << stats_.skeleton_voxels_initial << "\n";
        f << "skeleton_voxels_pruned=" << stats_.skeleton_voxels_pruned << "\n";
        f << "graph_node_count=" << stats_.graph_node_count << "\n";
        f << "graph_edge_count=" << stats_.graph_edge_count << "\n";
        f << "surface_node_count=" << stats_.surface_node_count << "\n";
        f << "interior_node_count=" << stats_.interior_node_count << "\n";
        f << "avg_degree=" << stats_.avg_degree << "\n";
        f << "final_vp=" << stats_.final_vp << "\n";
        f << "weights: w_angle=" << vp_weights_.w_angle
            << " w_length=" << vp_weights_.w_length
            << " w_location=" << vp_weights_.w_location
            << " w_direction=" << vp_weights_.w_direction << "\n";
        std::cout << "  Exported: " << fn << std::endl;
    }
}

// ============================================================================
// run() - execute full pipeline
// ============================================================================
PipelineStats TpmsPipeline::run(uint32_t steps_mask) {
    std::cout << "=== TPMS Pipeline ===" << std::endl;
    stats_.resolution = params_.resolution;
    stats_.bbox_half = params_.half_bbox;
    stats_.is_perturbed = params_.is_perturbed;

    auto t0 = std::chrono::high_resolution_clock::now();

    if (steps_mask == 0 || (steps_mask & 1))  step1_generateTpmsField();
    if (steps_mask == 0 || (steps_mask & 2))  step2_voxelize();
    if (steps_mask == 0 || (steps_mask & 4))  step3_extractLargestVoidRegion();
    if (steps_mask == 0 || (steps_mask & 8))  step4_skeletonize();
    if (steps_mask == 0 || (steps_mask & 16)) step5_pruneSkeleton();
    if (steps_mask == 0 || (steps_mask & 32)) step6_buildGraph();
    if (steps_mask == 0 || (steps_mask & 64)) step7_computeVP();
    if (steps_mask == 0 || (steps_mask & 128)) step8_exportData("tpms_output/default");

    auto t1 = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    std::cout << "=== Pipeline completed in " << dur.count() / 1000.0 << " s ===" << std::endl;

    return stats_;
}

// ============================================================================
// runTpmsComparison - run two cases and compare
// ============================================================================
void runTpmsComparison(int resolution, const string& output_dir) {
    random_device rd;
    mt19937 rng(rd());
    uniform_real_distribution<double> alpha_noise(0.9, 1.15);
    uniform_real_distribution<double> phase_noise(-0.25, 0.25);

    // ---- Case 1: Regular Gyroid ----
    std::cout << "\n========== CASE 1: Regular Gyroid ==========\n";
    {
        TpmsParams p;
        p.resolution = resolution;
        p.alpha_x = 1.0; p.alpha_y = 1.0; p.alpha_z = 1.0;
        p.phase_x = 0.0; p.phase_y = 0.0; p.phase_z = 0.0;
        p.is_perturbed = false;

        TpmsPipeline pipe;
        pipe.setParams(p);
        PipelineStats s = pipe.run();
        pipe.step8_exportData(output_dir + "/regular");

        std::cout << "\n--- Regular Gyroid Stats ---\n";
        std::cout << "  Graph nodes: " << s.graph_node_count << std::endl;
        std::cout << "  Graph edges: " << s.graph_edge_count << std::endl;
        std::cout << "  Surface nodes: " << s.surface_node_count << std::endl;
        std::cout << "  Avg degree: " << s.avg_degree << std::endl;
        std::cout << "  Final VP: " << s.final_vp << std::endl;
        std::cout << "  Skeleton voxels: " << s.skeleton_voxels_pruned << std::endl;
    }

    // ---- Case 2: Perturbed Gyroid ----
    std::cout << "\n========== CASE 2: Perturbed Gyroid ==========\n";
    {
        TpmsParams p;
        p.resolution = resolution;
        p.alpha_x = alpha_noise(rng); p.alpha_y = alpha_noise(rng); p.alpha_z = alpha_noise(rng);
        p.phase_x = phase_noise(rng); p.phase_y = phase_noise(rng); p.phase_z = phase_noise(rng);
        p.is_perturbed = true;

        TpmsPipeline pipe;
        pipe.setParams(p);
        PipelineStats s = pipe.run();
        pipe.step8_exportData(output_dir + "/perturbed");

        std::cout << "\n--- Perturbed Gyroid Stats ---\n";
        std::cout << "  alpha=(" << p.alpha_x << "," << p.alpha_y << "," << p.alpha_z << ")" << std::endl;
        std::cout << "  phase=(" << p.phase_x << "," << p.phase_y << "," << p.phase_z << ")" << std::endl;
        std::cout << "  Graph nodes: " << s.graph_node_count << std::endl;
        std::cout << "  Graph edges: " << s.graph_edge_count << std::endl;
        std::cout << "  Surface nodes: " << s.surface_node_count << std::endl;
        std::cout << "  Avg degree: " << s.avg_degree << std::endl;
        std::cout << "  Final VP: " << s.final_vp << std::endl;
        std::cout << "  Skeleton voxels: " << s.skeleton_voxels_pruned << std::endl;
    }

    std::cout << "\n========== Comparison Complete ==========\n";
}