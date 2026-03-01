#include "MorseComplex.h"
#include <Eigen/Dense>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>

using Vector3d = Eigen::Vector3d;

namespace {

    // ---------- grid helpers ----------
    inline int idx3(int x, int y, int z, int N) { return x * N * N + y * N + z; }

    inline bool in_bounds(int x, int y, int z, int N) {
        return (x >= 0 && x < N && y >= 0 && y < N && z >= 0 && z < N);
    }

    static const int k6[6][3] = {
        {+1,0,0},{-1,0,0},{0,+1,0},{0,-1,0},{0,0,+1},{0,0,-1}
    };

    static std::vector<int> neighbors26(int id, int N) {
        int x = id / (N * N);
        int r = id % (N * N);
        int y = r / N;
        int z = r % N;
        std::vector<int> out;
        out.reserve(26);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    if (dx == 0 && dy == 0 && dz == 0) continue;
                    int nx = x + dx, ny = y + dy, nz = z + dz;
                    if (in_bounds(nx, ny, nz, N)) out.push_back(idx3(nx, ny, nz, N));
                }
        return out;
    }

    static std::vector<int> neighbors6(int id, int N) {
        int x = id / (N * N);
        int r = id % (N * N);
        int y = r / N;
        int z = r % N;
        std::vector<int> out;
        out.reserve(6);
        for (auto& d : k6) {
            int nx = x + d[0], ny = y + d[1], nz = z + d[2];
            if (in_bounds(nx, ny, nz, N)) out.push_back(idx3(nx, ny, nz, N));
        }
        return out;
    }

    // 为了稳定处理 plateau/tie：加入极小扰动（按 id 打破平局）
    inline double f_tie(const Eigen::VectorXd& f, int id) {
        // 1e-16 级别扰动，一般不会改变真实数值关系，只用于避免 == 的歧义
        return f[id] + 1e-16 * (double)id;
    }

    static bool is_strict_max(int id, int N, const Eigen::VectorXd& f, double eps) {
        const double v0 = f_tie(f, id);
        for (int nb : neighbors26(id, N)) {
            if (f_tie(f, nb) >= v0 - eps) return false;
        }
        return true;
    }

    static bool is_strict_min(int id, int N, const Eigen::VectorXd& f, double eps) {
        const double v0 = f_tie(f, id);
        for (int nb : neighbors26(id, N)) {
            if (f_tie(f, nb) <= v0 + eps) return false;
        }
        return true;
    }

    // 计算 upper/lower link 的连通分量：
    // - link 节点来自 26 邻域
    // - 连通性用 6 邻接（与你文件 Is1SaddlePoint/Is2SaddlePoint 的 floodfill 风格一致）
    //static std::vector<std::vector<int>> connected_components_in_link(
    //    int id, int N, const Eigen::VectorXd& f, bool upper_link, double eps)
    //{
    //    const double v0 = f_tie(f, id);
    //    auto n26 = neighbors26(id, N);

    //    std::vector<char> mark(N * N * N, 0); // 0: out, 1: in_link(unvisited), 2: visited
    //    std::vector<int> link_nodes;
    //    link_nodes.reserve(n26.size());

    //    for (int nb : n26) {
    //        const double vnb = f_tie(f, nb);
    //        if (upper_link) {
    //            if (vnb > v0 + eps) { mark[nb] = 1; link_nodes.push_back(nb); }
    //        }
    //        else {
    //            if (vnb < v0 - eps) { mark[nb] = 1; link_nodes.push_back(nb); }
    //        }
    //    }

    //    std::vector<std::vector<int>> comps;
    //    std::queue<int> q;

    //    for (int seed : link_nodes) {
    //        if (mark[seed] != 1) continue;
    //        comps.emplace_back();
    //        mark[seed] = 2;
    //        q.push(seed);
    //        while (!q.empty()) {
    //            int cur = q.front(); q.pop();
    //            comps.back().push_back(cur);
    //            for (int nb6 : neighbors6(cur, N)) {
    //                if (mark[nb6] == 1) {
    //                    mark[nb6] = 2;
    //                    q.push(nb6);
    //                }
    //            }
    //        }
    //    }
    //    return comps;
    //}

    static std::vector<std::vector<int>> connected_components_in_link(
        int id, int N, const Eigen::VectorXd& f,
        bool upper_link, double eps)
    {
        const double v0 = f_tie(f, id);
        auto n26 = neighbors26(id, N);

        std::unordered_set<int> link_set;
        link_set.reserve(26);

        for (int nb : n26)
        {
            double vnb = f_tie(f, nb);
            if (upper_link)
            {
                if (vnb > v0 + eps)
                    link_set.insert(nb);
            }
            else
            {
                if (vnb < v0 - eps)
                    link_set.insert(nb);
            }
        }

        std::vector<std::vector<int>> comps;
        std::unordered_set<int> visited;
        visited.reserve(26);

        for (int seed : link_set)
        {
            if (visited.count(seed)) continue;

            comps.emplace_back();
            std::queue<int> q;
            q.push(seed);
            visited.insert(seed);

            while (!q.empty())
            {
                int cur = q.front(); q.pop();
                comps.back().push_back(cur);

                for (int nb : neighbors6(cur, N))
                {
                    if (link_set.count(nb) && !visited.count(nb))
                    {
                        visited.insert(nb);
                        q.push(nb);
                    }
                }
            }
        }

        return comps;
    }



    static bool is_2_saddle(int id, int N, const Eigen::VectorXd& f, double eps) {
        auto comps = connected_components_in_link(id, N, f, /*upper_link=*/true, eps);
        return comps.size() >= 2;
    }
    static bool is_1_saddle(int id, int N, const Eigen::VectorXd& f, double eps) {
        auto comps = connected_components_in_link(id, N, f, /*upper_link=*/false, eps);
        return comps.size() >= 2;
    }

    // 离散梯度追踪：从 start 沿最陡上升/下降走到"无法继续"的点
    static int follow_discrete_gradient(int start, int N, const Eigen::VectorXd& f,
        bool ascent, double eps, int max_steps = 100000)
    {
        int cur = start;
        for (int step = 0; step < max_steps; ++step) {
            int best_id = cur;
            double best = f_tie(f, cur);

            for (int nb : neighbors26(cur, N)) {
                double v = f_tie(f, nb);
                if (ascent) {
                    if (v > best + eps) { best = v; best_id = nb; }
                }
                else {
                    if (v < best - eps) { best = v; best_id = nb; }
                }
            }
            if (best_id == cur) break;
            cur = best_id;
        }
        return cur;
    }

    // 最近 kernel（简单 O(K)；若 kernel 很多再换 kd-tree）
    static int nearest_kernel(const std::vector<Vector3d>& kernels, const Vector3d& p) {
        int best = -1;
        double best_d2 = std::numeric_limits<double>::infinity();
        for (int i = 0; i < (int)kernels.size(); ++i) {
            double d2 = (kernels[i] - p).squaredNorm();
            if (d2 < best_d2) { best_d2 = d2; best = i; }
        }
        return best;
    }

    struct PairHash {
        std::size_t operator()(const std::pair<int, int>& p) const noexcept {
            // 64-bit hash combine
            return (std::size_t)((uint64_t)(uint32_t)p.first << 32) ^ (uint32_t)p.second;
        }
    };

} // namespace


std::vector<std::pair<int, int>>  MorseComplex::compare_msc(std::vector<Vector3d> Kernels,
    Eigen::VectorXd SDF_gaussian,
    int res,
    int grid_num)
{
    // -------- 0) 校验 N --------
    const int M = (int)SDF_gaussian.size();
    const int N = res;

     std::cerr << "[compare_msc] size mismatch: res=" << res
            << ", res^3=" << (long long)res * res * res
            << ", SDF_gaussian.size()=" << M
            << ", grid_num=" << grid_num << "\n";

    // -------- 1) 坐标系：[-0.5,0.5]^3 --------
    const double xmin = -0.5;
    const double xmax = +0.5;
    const double h = (N > 1) ? (xmax - xmin) / (double)(N - 1) : 0.0;

    auto grid_pos = [&](int id) -> Vector3d {
        int x = id / (N * N);
        int r = id % (N * N);
        int y = r / N;
        int z = r % N;
        return Vector3d(xmin + x * h, xmin + y * h, xmin + z * h);
        };

    // -------- 2) 扫描临界点（max/min/1-saddle/2-saddle）--------
    const double eps = 1e-14; // tie-break 后 eps 可更小
    std::vector<int> maxima, minima, saddle1, saddle2;
    maxima.reserve(M / 200);
    minima.reserve(M / 200);
    saddle1.reserve(M / 200);
    saddle2.reserve(M / 200);

    for (int id = 0; id < M; ++id) {
        int x = id / (N * N);
        int r = id % (N * N);
        int y = r / N;
        int z = r % N;

        // 边界点很容易出现伪临界点：先跳过最外层
        if (x == 0 || y == 0 || z == 0 || x == N - 1 || y == N - 1 || z == N - 1) continue;

        if (is_strict_max(id, N, SDF_gaussian, eps)) maxima.push_back(id);
        if (is_strict_min(id, N, SDF_gaussian, eps)) minima.push_back(id);

        if (is_2_saddle(id, N, SDF_gaussian, eps)) saddle2.push_back(id);
        if (is_1_saddle(id, N, SDF_gaussian, eps)) saddle1.push_back(id);
    }

    std::cout << "[MSC] grid N=" << N << " range=[-0.5,0.5]^3\n";
    std::cout << "[MSC] maxima=" << maxima.size()
        << ", minima=" << minima.size()
        << ", 1-saddle=" << saddle1.size()
        << ", 2-saddle=" << saddle2.size() << "\n";

    //if (Kernels.empty()) {
    //    std::cerr << "[MSC] Kernels is empty.\n";
    //    return;
    //}

    // -------- 3) maximum -> kernel 映射（对应你文件里 maximum 的 sSourceID）--------
    std::unordered_map<int, int> max_to_kernel;
    max_to_kernel.reserve(maxima.size() * 2);

    for (int mid : maxima) {
        int kid = nearest_kernel(Kernels, grid_pos(mid));
        max_to_kernel[mid] = kid;
    }

    // -------- 4) 用 2-saddle 构建 kernel 图（去重：保留 saddle 值更大者，类似 reduce_pairs for 2-saddle）--------
    // edge_best[(k1,k2)] = best saddle value
    std::unordered_map<std::pair<int, int>, double, PairHash> edge_best;
    int raw_pairs = 0;

    for (int sid : saddle2) {
        // upper link 分量
        auto comps = connected_components_in_link(sid, N, SDF_gaussian, /*upper_link=*/true, eps);
        if (comps.size() < 2) continue;

        // 对每个分量：取分量内 f 最大点作为 seed，沿 ascent 追踪到 max
        std::vector<int> reached_max;
        reached_max.reserve(comps.size());

        for (auto& comp : comps) {
            int seed = comp[0];
            for (int v : comp) {
                if (f_tie(SDF_gaussian, v) > f_tie(SDF_gaussian, seed)) seed = v;
            }
            int end = follow_discrete_gradient(seed, N, SDF_gaussian, /*ascent=*/true, eps);
            if (is_strict_max(end, N, SDF_gaussian, eps)) reached_max.push_back(end);
        }

        // 去重
        std::sort(reached_max.begin(), reached_max.end());
        reached_max.erase(std::unique(reached_max.begin(), reached_max.end()), reached_max.end());
        if (reached_max.size() < 2) continue;

        // 若出现 >2 个上升分量（离散误差可能发生）：选"max 值最高"的两个
        std::sort(reached_max.begin(), reached_max.end(), [&](int a, int b) {
            return f_tie(SDF_gaussian, a) > f_tie(SDF_gaussian, b);
            });
        int m1 = reached_max[0];
        int m2 = reached_max[1];

        int k1 = max_to_kernel.count(m1) ? max_to_kernel[m1] : -1;
        int k2 = max_to_kernel.count(m2) ? max_to_kernel[m2] : -1;
        if (k1 < 0 || k2 < 0 || k1 == k2) continue;

        if (k1 > k2) std::swap(k1, k2);
        ++raw_pairs;

        const double saddle_value = SDF_gaussian[sid]; // 2-saddle 的标量值
        auto key = std::make_pair(k1, k2);

        auto it = edge_best.find(key);
        if (it == edge_best.end()) {
            edge_best[key] = saddle_value;
        }
        else {
            // 2-saddle：保留更"高"的 saddle（对应你 reduce_pairs 对 2-saddle 的逻辑）
            if (saddle_value > it->second) it->second = saddle_value;
        }
    }

    std::cout << "[MSC] raw 2-saddle kernel pairs (before dedup): " << raw_pairs << "\n";
    std::cout << "[MSC] unique kernel edges (after dedup): " << edge_best.size() << "\n";

    std::vector<std::pair<int, int>> pairs;
    pairs.reserve(edge_best.size());  // 预分配内存以提高效率

    for (const auto& kv : edge_best) {
        pairs.push_back(kv.first);  // kv.first 就是 pair<int, int>
    }

	return pairs;
    // -------- 5) 生成邻接表 & 连通分量 --------
    //std::vector<std::vector<int>> adj(Kernels.size());
    //adj.assign(Kernels.size(), {});

    //for (auto& kv : edge_best) {
    //    int a = kv.first.first;
    //    int b = kv.first.second;
    //    adj[a].push_back(b);
    //    adj[b].push_back(a);
    //}
    //for (auto& v : adj) {
    //    std::sort(v.begin(), v.end());
    //    v.erase(std::unique(v.begin(), v.end()), v.end());
    //}

    //// BFS components
    //std::vector<char> vis(Kernels.size(), 0);
    //int comp_cnt = 0;
    //std::vector<int> comp_size;

    //for (int i = 0; i < (int)Kernels.size(); ++i) {
    //    if (vis[i]) continue;
    //    ++comp_cnt;
    //    int sz = 0;
    //    std::queue<int> q;
    //    q.push(i);
    //    vis[i] = 1;
    //    while (!q.empty()) {
    //        int u = q.front(); q.pop();
    //        ++sz;
    //        for (int v : adj[u]) {
    //            if (!vis[v]) { vis[v] = 1; q.push(v); }
    //        }
    //    }
    //    comp_size.push_back(sz);
    //}

    //std::sort(comp_size.begin(), comp_size.end(), std::greater<int>());

    //std::cout << "[MSC] kernel graph components = " << comp_cnt << "\n";
    //if (!comp_size.empty()) {
    //    std::cout << "[MSC] component sizes (desc): ";
    //    for (int i = 0; i < (int)std::min<size_t>(comp_size.size(), 10); ++i)
    //        std::cout << comp_size[i] << (i + 1 < (int)std::min<size_t>(comp_size.size(), 10) ? ", " : "");
    //    std::cout << "\n";
    //}

    // -------- 6) 输出 edge list（可用来 compare/可视化）--------
    // 这里直接打印前若干条；你也可以改成写文件（CSV/JSON）
    //int shown = 0;
    //std::cout << "[MSC] sample edges (k1,k2,saddle_value):\n";
    //for (auto& kv : edge_best) {
    //    if (shown >= 20) break;
    //    std::cout << "  (" << kv.first.first << ", " << kv.first.second
    //        << ", " << kv.second << ")\n";
    //    ++shown;
    //}
}