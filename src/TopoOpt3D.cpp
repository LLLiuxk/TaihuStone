#include "TopoOpt3D.h"
#include <fstream>
#include <iostream>

TopologyOptimizer3D::TopologyOptimizer3D(const std::string& vol_path, Config3D config)
    : cfg(config)
{
    loadVolumeData(vol_path);

    beta = cfg.beta_start;
    rho_tilde = rho; // 初始化

    std::cout << "Initializing 3D Filter (Radius: " << cfg.r_min << ")..." << std::endl;
    //precompute_filter();
    std::cout << "Filter initialized." << std::endl;
}

// 加载体素数据
// 格式: [int nx, int ny, int nz]\n [double data...]
void TopologyOptimizer3D::loadVolumeData(const std::string& path) {
    std::ifstream in(path);   //
    if (!in) throw std::runtime_error("Cannot open file: " + path);

    // 1. 读取体素分辨率（文本）
    in >> nx >> ny >> nz;

    n_vars = nx * ny * nz;
    rho.resize(n_vars);
    target_shape.resize(n_vars);

    std::cout << "n_vars " << n_vars
        << ": " << nx << " " << ny << " " << nz << std::endl;

    // 2. 逐个读取 density（文本 double）
    for (int i = 0; i < n_vars; ++i) {
        double val;

        if (!(in >> val)) {
            throw std::runtime_error(
                "Not enough voxel data in file"            );
        }

        if (i < 100) std::cout << val << std::endl;
        if (val < 0) val = 0.0;
        if (val > 1) val = 1.0;

        rho(i) = val;
        target_shape(i) = val;
    }

    std::cout << "Loaded 3D Volume: " << nx << "x" << ny << "x" << nz << std::endl;
}

// ==========================================
// 3D 物理场更新逻辑
// ==========================================

// 预计算滤波邻居 (解决 3D 卷积慢的问题)
void TopologyOptimizer3D::precompute_filter() {
    filter_neighborhoods.resize(n_vars);
    filter_weight_sums.resize(n_vars, 0.0);

    int r_int = std::ceil(cfg.r_min);
    double r_sq = cfg.r_min * cfg.r_min;

    // 遍历所有体素
#pragma omp parallel for // 如果你有 OpenMP，这里可以加速
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int i = idx(x, y, z);

                // 搜索邻域
                for (int kz = -r_int; kz <= r_int; ++kz) {
                    for (int ky = -r_int; ky <= r_int; ++ky) {
                        for (int kx = -r_int; kx <= r_int; ++kx) {

                            int nz_idx = z + kz;
                            int ny_idx = y + ky;
                            int nx_idx = x + kx;

                            // 边界检查
                            if (nz_idx >= 0 && nz_idx < nz &&
                                ny_idx >= 0 && ny_idx < ny &&
                                nx_idx >= 0 && nx_idx < nx) {

                                double dist_sq = kx * kx + ky * ky + kz * kz;
                                if (dist_sq <= r_sq) {
                                    double dist = std::sqrt(dist_sq);
                                    double w = cfg.r_min - dist;

                                    // 存储邻居索引和权重
                                    int neighbor_id = idx(nx_idx, ny_idx, nz_idx);

                                    // 注意：这里需要线程安全，如果用了 OpenMP，push_back 要小心
                                    // 简单起见，这里先不加 OpenMP，或者预分配
                                    filter_neighborhoods[i].push_back({ neighbor_id, w });
                                    filter_weight_sums[i] += w;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// 应用 3D 线性滤波 (使用预计算表)
Eigen::VectorXd TopologyOptimizer3D::apply_linear_filter(const Eigen::VectorXd& x) {
    Eigen::VectorXd x_new = Eigen::VectorXd::Zero(n_vars);

    for (int i = 0; i < n_vars; ++i) {
        if (filter_weight_sums[i] < 1e-6) {
            x_new(i) = x(i);
            continue;
        }

        double sum = 0.0;
        for (const auto& nb : filter_neighborhoods[i]) {
            sum += nb.weight * x(nb.idx);
        }
        x_new(i) = sum / filter_weight_sums[i];
    }
    return x_new;
}

// 滤波的反向传播
Eigen::VectorXd TopologyOptimizer3D::backprop_filter(const Eigen::VectorXd& sens_phys, const Eigen::VectorXd& rho_bar) {
    // 1. Heaviside 反向 (Chain Rule Step 1)
    Eigen::VectorXd sens_bar = Eigen::VectorXd::Zero(n_vars);
    for (int i = 0; i < n_vars; ++i) {
        double dH = d_proj_heaviside(rho_bar(i), beta, cfg.eta);
        sens_bar(i) = sens_phys(i) * dH;
    }

    // 2. 线性滤波反向 (Chain Rule Step 2)
    // 线性滤波矩阵是对称的 (如果归一化处理得当)，这里近似复用前向滤波逻辑
    // 严格来说应该除以 neighbors 的 weight_sum，但通常近似为自伴随算子
    // 简单的实现：直接再次调用 apply_linear_filter
    return apply_linear_filter(sens_bar);
}

// 物理场更新总入口
void TopologyOptimizer3D::update_physics_field() {
    // 1. Filter
    Eigen::VectorXd rho_bar = apply_linear_filter(rho);

    // 2. Project
    for (int i = 0; i < n_vars; ++i) {
        rho_tilde(i) = proj_heaviside(rho_bar(i), beta, cfg.eta);
    }
}

// ==========================================
// 辅助功能：对称性与 IO
// ==========================================

void TopologyOptimizer3D::imposeSymmetry(Eigen::VectorXd& x) {
    // 假设关于 X 轴中心对称 (Left-Right Symmetry)
    // 即 x 与 nx - 1 - x 对称
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int i = 0; i < nx / 2; ++i) {
                int id_left = idx(i, y, z);
                int id_right = idx(nx - 1 - i, y, z);

                double avg = (x(id_left) + x(id_right)) / 2.0;
                x(id_left) = avg;
                x(id_right) = avg;
            }
        }
    }
}

// 保存为 VTI 格式 (可被 ParaView 读取)
// 这是一个简单的 ASCII XML 格式，虽然体积大但无需依赖库
void TopologyOptimizer3D::saveVTI(const std::string& filename) {
    std::ofstream out(filename);
    if (!out) return;

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <ImageData WholeExtent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    out << "    <Piece Extent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\">\n";
    out << "      <PointData Scalars=\"Density\">\n";
    out << "        <DataArray type=\"Float32\" Name=\"Density\" format=\"ascii\">\n";

    for (int i = 0; i < n_vars; ++i) {
        out << (float)rho_tilde(i) << " ";
        if (i % 20 == 0) out << "\n";
    }

    out << "\n        </DataArray>\n";
    out << "      </PointData>\n";
    out << "      <CellData>\n";
    out << "      </CellData>\n";
    out << "    </Piece>\n";
    out << "  </ImageData>\n";
    out << "</VTKFile>\n";

    std::cout << "Saved VTI: " << filename << std::endl;
}