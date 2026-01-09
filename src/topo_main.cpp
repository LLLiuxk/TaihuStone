#include "TopoOpt.h"

int main(int argc, char** argv) {
    //initTShape();
    Config cfg;
    cfg.alpha_crit_deg = 45.0; // 45度悬垂角
    cfg.max_iter = 5000;         // 演示用
    cfg.r_min = 8.0;
    cfg.max_penalty = 3000.0;
    cfg.beta_max = 24.0;
    cfg.vol_frac = 0.5;

    std::string filename = "topology.jpg";
    TopologyOptimizer opt(filename, cfg);
    std::cout << "Starting Optimization..." << std::endl;
    opt.solve();
    std::cout << "Done. Result saved to output_optimized.png" << std::endl;

    return 0;
}