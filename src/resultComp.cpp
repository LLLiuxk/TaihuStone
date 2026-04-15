#include "resultComp.h"
#include <cmath>
#include <stack>
#include <algorithm>
#include <iostream>

// 如果编译器支持 OpenMP，则包含头文件
#ifdef _OPENMP
#include <omp.h>
#endif

GaussianFieldMSC::GaussianFieldMSC(int nx, int ny, int nz, Vec3 min_b, Vec3 max_b)
    : nx_(nx), ny_(ny), nz_(nz), min_b_(min_b), max_b_(max_b) {
    // 计算物理步长
    dx_ = (nx_ > 1) ? (max_b_.x - min_b_.x) / (nx_ - 1) : 0;
    dy_ = (ny_ > 1) ? (max_b_.y - min_b_.y) / (ny_ - 1) : 0;
    dz_ = (nz_ > 1) ? (max_b_.z - min_b_.z) / (nz_ - 1) : 0;
}

std::vector<std::vector<int>> GaussianFieldMSC::ComputeConnectivity(
    const Eigen::VectorXd& field_data,
    const std::vector<Vec3>& original_points)
{
    if (field_data.size() != nx_ * ny_ * nz_) {
        std::cerr << "Field size mismatch.\n";
        return {};
    }

    const int total_size = nx_ * ny_ * nz_;
    const double eps = 1e-12;

    // -------------------------------
    // 1. 找 maxima
    // -------------------------------
    std::vector<char> is_maxima(total_size, 0);

#pragma omp parallel for
    for (int i = 0; i < total_size; ++i) {
        int x, y, z;
        getCoord(i, x, y, z);
        if (x <= 0 || x >= nx_ - 1 ||
            y <= 0 || y >= ny_ - 1 ||
            z <= 0 || z >= nz_ - 1)
            continue;

        if (CheckIsMaxima(field_data, x, y, z))
            is_maxima[i] = 1;
    }

    // -------------------------------
    // 2. 找 2-saddles
    // -------------------------------
    std::vector<int> saddles;

#pragma omp parallel
    {
        std::vector<int> local_saddles;

#pragma omp for nowait
        for (int i = 0; i < total_size; ++i) {
            int x, y, z;
            getCoord(i, x, y, z);
            if (x <= 0 || x >= nx_ - 1 ||
                y <= 0 || y >= ny_ - 1 ||
                z <= 0 || z >= nz_ - 1)
                continue;

            if (CheckIs2Saddle(field_data, x, y, z))
                local_saddles.push_back(i);
        }

#pragma omp critical
        saddles.insert(saddles.end(),
            local_saddles.begin(),
            local_saddles.end());
    }

    // -------------------------------
    // 3. 构建边（折中策略）
    // -------------------------------

    // 用 map 存 pair -> best saddle value
    std::unordered_map<long long, double> edge_best;

    for (int saddle_idx : saddles)
    {
        auto reached = TraceAscent(field_data,
            reinterpret_cast<const std::vector<bool>&>(is_maxima),
            saddle_idx);

        if (reached.size() < 2)
            continue;

        // 转为 vector
        std::vector<int> maxima_list(reached.begin(), reached.end());

        // 按 field 值排序
        std::sort(maxima_list.begin(), maxima_list.end(),
            [&](int a, int b)
            {
                return field_data[a] > field_data[b];
            });

        // 只取前两个
        int m1 = maxima_list[0];
        int m2 = maxima_list[1];

        double saddle_value = field_data[saddle_idx];

        // ---- Persistence 计算 ----
        double persistence =
            std::min(field_data[m1], field_data[m2])
            - saddle_value;

        // 阈值 tau
        double f_min = field_data.minCoeff();
        double f_max = field_data.maxCoeff();
        double global_range = f_max - f_min;

        double persistence_threshold = 0.9 * global_range;
        double tau = persistence_threshold;

        // 过滤弱连接
        if (persistence < tau)
            continue;

        int p1 = MapGridToPoint(m1, original_points);
        int p2 = MapGridToPoint(m2, original_points);

        if (p1 < 0 || p2 < 0 || p1 == p2)
            continue;

        if (p1 > p2) std::swap(p1, p2);

        long long key = ((long long)p1 << 32) | p2;

        auto it = edge_best.find(key);
        if (it == edge_best.end())
            edge_best[key] = saddle_value;
        else
            if (saddle_value > it->second)
                it->second = saddle_value;
    }

    // -------------------------------
    // 4. 转 adjacency
    // -------------------------------
    std::vector<std::vector<int>> adjacency(original_points.size());

    for (auto& kv : edge_best)
    {
        long long key = kv.first;
        int u = (int)(key >> 32);
        int v = (int)(key & 0xffffffff);

        adjacency[u].push_back(v);
        adjacency[v].push_back(u);
    }

    return adjacency;
}


//std::vector<std::vector<int>> GaussianFieldMSC::ComputeConnectivity(
//    const Eigen::VectorXd& field_data,
//    const std::vector<Vec3>& original_points)
//{
//    // 0. 数据校验
//    if (field_data.size() != nx_ * ny_ * nz_) {
//        std::cerr << "[GaussianFieldMSC] Error: Field data size mismatch! "
//            << "Expected " << nx_ * ny_ * nz_ << ", got " << field_data.size() << std::endl;
//        return {};
//    }
//
//    // 1. 识别临界点 (Maxima 和 2-Saddles)
//    std::vector<int> saddles;
//    // 使用 std::vector<char> 代替 bool 以避免某些 vector<bool> 的并发写入问题，或者仅作为只读
//    // 这里我们先串行初始化
//    std::vector<bool> is_maxima(field_data.size(), false);
//
//    // 并行寻找临界点
//    // 注意：vector<bool> 不是线程安全的写入，所以我们这里分两步，或者使用 int 数组
//    // 为了简单起见，这里先计算 Maxima 标记，再计算 Saddle
//
//    // Step 1.1: 标记 Maxima
//#pragma omp parallel for
//    for (int z = 1; z < nz_ - 1; ++z) {
//        for (int y = 1; y < ny_ - 1; ++y) {
//            for (int x = 1; x < nx_ - 1; ++x) {
//                if (CheckIsMaxima(field_data, x, y, z)) {
//                    int idx = getIdx(x, y, z);
//                    // 即使是 vector<bool>，不同索引的并发写入通常是不安全的(bit操作)
//                    // 但在这里我们为了严谨，应该避免并发写 vector<bool>。
//                    // 考虑到 Maxima 稀疏，我们这里暂不并行写入，或者改用 vector<char>。
//                    // 修正方案：此处不直接写全局 vector<bool>，或者取消并行。
//                    // 鉴于性能瓶颈通常在 Saddle 判断，Maxima 判断很快，这里取消并行写入，
//                    // 或者使用 critical (虽然慢)。
//                    // 最好的方式：
//                }
//            }
//        }
//    }
//
//    // 重新实现并行安全的 Maxima 检测
//    std::vector<char> is_maxima_char(field_data.size(), 0);
//#pragma omp parallel for
//    for (int i = 0; i < nx_ * ny_ * nz_; ++i) {
//        int x, y, z;
//        getCoord(i, x, y, z);
//        // 排除边界
//        if (x > 0 && x < nx_ - 1 && y > 0 && y < ny_ - 1 && z > 0 && z < nz_ - 1) {
//            if (CheckIsMaxima(field_data, x, y, z)) {
//                is_maxima_char[i] = 1;
//            }
//        }
//    }
//    // 转回 bool vector 供后续使用 (其实直接用 char 也可以，为了匹配接口转一下)
//    for (size_t i = 0; i < is_maxima.size(); ++i) is_maxima[i] = (is_maxima_char[i] == 1);
//
//    // Step 1.2: 寻找 Saddles
//#pragma omp parallel for
//    for (int z = 1; z < nz_ - 1; ++z) {
//        for (int y = 1; y < ny_ - 1; ++y) {
//            for (int x = 1; x < nx_ - 1; ++x) {
//                if (CheckIs2Saddle(field_data, x, y, z)) {
//                    int idx = getIdx(x, y, z);
//#pragma omp critical
//                    {
//                        saddles.push_back(idx);
//                    }
//                }
//            }
//        }
//    }
//
//    // std::cout << "[GaussianFieldMSC] Found " << saddles.size() << " 2-Saddles." << std::endl;
//
//    // 2. 追踪脊线 (Separatrix) 并建立连接
//    std::vector<std::vector<int>> adjacency(original_points.size());
//
//    // 这部分也可以并行，但需要注意 adjacency 的写入安全
//    // 考虑到 saddles 数量通常不多，串行即可，或者使用局部归约
//    for (int saddle_idx : saddles) {
//        // 从鞍点出发，沿着梯度上升，找到所有能到达的极大值点
//        std::set<int> connected_grid_maxima = TraceAscent(field_data, is_maxima, saddle_idx);
//
//        // 如果一个鞍点连接了 >= 2 个极大值，说明这两个极大值在拓扑上相邻
//        if (connected_grid_maxima.size() >= 2) {
//            // 将网格上的极大值点索引，转换为用户输入的 original_points 的索引
//            std::vector<int> mapped_indices;
//            for (int grid_max_idx : connected_grid_maxima) {
//                int p_idx = MapGridToPoint(grid_max_idx, original_points);
//                if (p_idx != -1) mapped_indices.push_back(p_idx);
//            }
//
//            // 在这些点之间建立两两连接
//            for (size_t i = 0; i < mapped_indices.size(); ++i) {
//                for (size_t j = i + 1; j < mapped_indices.size(); ++j) {
//                    int u = mapped_indices[i];
//                    int v = mapped_indices[j];
//                    if (u == v) continue;
//
//                    // 添加无向边 (防止重复)
//                    bool exists = false;
//                    for (int neighbor : adjacency[u]) if (neighbor == v) exists = true;
//                    if (!exists) adjacency[u].push_back(v);
//
//                    exists = false;
//                    for (int neighbor : adjacency[v]) if (neighbor == u) exists = true;
//                    if (!exists) adjacency[v].push_back(u);
//                }
//            }
//        }
//    }
//
//    return adjacency;
//}

// ---------------- Private Helper Functions ----------------

inline int GaussianFieldMSC::getIdx(int x, int y, int z) const {
    return z * ny_ * nx_ + y * nx_ + x;
}

inline void GaussianFieldMSC::getCoord(int idx, int& x, int& y, int& z) const {
    x = idx % nx_;
    y = (idx / nx_) % ny_;
    z = idx / (nx_ * ny_);
}

bool GaussianFieldMSC::CheckIsMaxima(const Eigen::VectorXd& field, int x, int y, int z) {
    double val = field[getIdx(x, y, z)];
    // 检查 26 邻域
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) continue;
                if (field[getIdx(x + dx, y + dy, z + dz)] > val) return false;
            }
        }
    }
    return true;
}

bool GaussianFieldMSC::CheckIs2Saddle(const Eigen::VectorXd& field, int x, int y, int z) {
    int center_idx = getIdx(x, y, z);
    double center_val = field[center_idx];

    std::vector<int> upper_neighbors;
    // 1. 收集所有比中心点值大的邻居 (Upper Star)
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) continue;
                int n_idx = getIdx(x + dx, y + dy, z + dz);
                if (field[n_idx] > center_val) {
                    upper_neighbors.push_back(n_idx);
                }
            }
        }
    }

    if (upper_neighbors.empty()) return false; // 极大值点，不是鞍点

    // 2. 计算连通分量 (使用 BFS)
    int components = 0;
    std::set<int> visited; // 记录已访问的 upper_neighbors

    for (int start_node : upper_neighbors) {
        if (visited.count(start_node)) continue;

        components++;
        if (components > 1) return true; // 只要分量 > 1，即为鞍点

        // BFS 遍历当前分量
        std::stack<int> s;
        s.push(start_node);
        visited.insert(start_node);

        while (!s.empty()) {
            int curr = s.top(); s.pop();
            int cx, cy, cz;
            getCoord(curr, cx, cy, cz);

            // 检查 curr 的邻居是否也在 upper_neighbors 中且空间相邻
            // 这里使用 6-邻域连通性 (更严格，有助于分离分量)
            // 如果使用 26-邻域，某些对角线连接可能会合并本应分开的分量
            for (int dz = -1; dz <= 1; ++dz) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        // 仅允许直接相邻 (Manhattan distance = 1)
                        if (std::abs(dx) + std::abs(dy) + std::abs(dz) != 1) continue;

                        int nx = cx + dx, ny = cy + dy, nz = cz + dz;
                        // 边界检查
                        if (nx < 0 || nx >= nx_ || ny < 0 || ny >= ny_ || nz < 0 || nz >= nz_) continue;

                        int n_idx = getIdx(nx, ny, nz);

                        // 必须是“高值邻居”集合中的一员
                        bool is_in_upper = false;
                        for (int un : upper_neighbors) if (un == n_idx) { is_in_upper = true; break; }

                        if (is_in_upper && visited.find(n_idx) == visited.end()) {
                            visited.insert(n_idx);
                            s.push(n_idx);
                        }
                    }
                }
            }
        }
    }
    return components >= 2;
}

std::set<int> GaussianFieldMSC::TraceAscent(const Eigen::VectorXd& field, const std::vector<bool>& is_maxima, int start_idx) {
    std::set<int> reached_maxima;

    int x, y, z;
    getCoord(start_idx, x, y, z);
    double center_val = field[start_idx];

    // 鞍点周围可能有多个上升方向，分别追踪
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) continue;

                int n_idx = getIdx(x + dx, y + dy, z + dz);

                // 如果邻居比鞍点高，则沿着这个方向爬山
                if (field[n_idx] > center_val) {
                    int curr = n_idx;
                    // 贪婪爬山
                    while (!is_maxima[curr]) {
                        int next_best = -1;
                        double max_v = field[curr];

                        int cx, cy, cz;
                        getCoord(curr, cx, cy, cz);

                        // 在 26-邻域中找最大值
                        for (int k = -1; k <= 1; ++k) {
                            for (int j = -1; j <= 1; ++j) {
                                for (int i = -1; i <= 1; ++i) {
                                    if (i == 0 && j == 0 && k == 0) continue;
                                    // 边界检查
                                    if (cx + i <= 0 || cx + i >= nx_ - 1 || cy + j <= 0 || cy + j >= ny_ - 1 || cz + k <= 0 || cz + k >= nz_ - 1) continue;

                                    int neighbor = getIdx(cx + i, cy + j, cz + k);
                                    if (field[neighbor] > max_v) {
                                        max_v = field[neighbor];
                                        next_best = neighbor;
                                    }
                                }
                            }
                        }

                        if (next_best != -1) {
                            curr = next_best;
                        }
                        else {
                            // 到达平坦区或局部极值，停止
                            break;
                        }
                    }
                    reached_maxima.insert(curr);
                }
            }
        }
    }
    return reached_maxima;
}

int GaussianFieldMSC::MapGridToPoint(int grid_idx, const std::vector<Vec3>& points) {
    int x, y, z;
    getCoord(grid_idx, x, y, z);

    // 计算网格点的物理坐标
    Vec3 pos;
    pos.x = min_b_.x + x * dx_;
    pos.y = min_b_.y + y * dy_;
    pos.z = min_b_.z + z * dz_;

    int best_idx = -1;
    double min_dist_sq = 1e20;

    for (size_t i = 0; i < points.size(); ++i) {
        double d = pos.distSq(points[i]);
        if (d < min_dist_sq) {
            min_dist_sq = d;
            best_idx = i;
        }
    }
    return best_idx;
}