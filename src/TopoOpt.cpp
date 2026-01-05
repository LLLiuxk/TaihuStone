#include "GeometryCorrector.h"

double Match_weight = 200.0;

double R_min = 8.0;

double Vol_weight = 0.1;

double Overhang_weight = 20.0;

double Vol_frac = 0.45;

double Max_beta = 24.0;

int main() {
    // 1. 初始化参数
    double init_learning_rate = 0.05;
    int logical_w = 40; // 网格宽度
    int logical_h = 20; // 网格高度
    double r_min = R_min; // 滤波半径
    int padding = 1;    // 设置填充
    //double vol_frac = 0.5; // 目标体积比 40%

    GeometryCorrector solver(logical_w, logical_h, padding, r_min);
    //GeometryCorrector solver(nx, ny, r_min);
    std::vector<double>& x = solver.getDesignVariables();
    //std::vector<double> target_mask = solver.initTShape(nx, ny, x);
    //solver.initializeTestShape();
    solver.initTShape();
    solver.printTopology(); // 打印当前拓扑
    std::cout << "Starting Topology Optimization with Overhang Constraints..." << std::endl;

    // 2. 优化主循环
    int max_iter = 20000;
    double learning_rate = 0.05; // 学习率

    std::vector<double> x_old = solver.getDesignVariables();
    double change_tolerance = 0.002; // 收敛容差 (0.5% 的变化)

    for (int iter = 0; iter < max_iter; ++iter) {
        // A. 物理场计算
        solver.applyFilter();
        solver.applyProjection();

        // B. 目标函数 (非对称保真)
        double obj = solver.computeObjectiveAndGradient();

        // C. 悬垂约束 (强力驱动)
        // 使用较大的权重来驱动生长
        double base_overhang = Overhang_weight;
        double overhang_weight = base_overhang + (iter / 50.0) * 5.0;
        if (overhang_weight > 100.0) overhang_weight = 100.0;

        //double overhang_weight = Overhang_weight;
        solver.applyOverhangConstraint(overhang_weight);

        // D. 链式法则
        solver.applyChainRule();

        // E. --- 关键：重置 Padding 的梯度 ---
        // 保证边界条件不被破坏
        solver.resetBoundarySensitivities();

        std::vector<double>& x = solver.getDesignVariables();
        const std::vector<double>& final_sens = solver.getSensitivities();
        x_old = x;
        // F. 更新设计变量 (梯度下降)
        //gradientDescentUpdate(x, final_sens, learning_rate);
        
        // F. 更新设计变量 (使用当前的 current_lr)
                // --- 动态调整学习率 ---
        // 策略：随着迭代进行，指数级降低学习率
        // 效果：前期 (iter < 500) 步子大，快速变形；后期 (iter > 1000) 步子极其细碎，用于微调
        double current_lr = init_learning_rate * std::exp(-0.002 * iter);

        // 设置一个下限，防止完全不动
        if (current_lr < 0.002) current_lr = 0.002;

        gradientDescentUpdate(x, final_sens, current_lr); // <--- 注意这里传入 current_lr
        solver.imposeSymmetry();
        
        // G. Beta 延拓
        if (iter > 0 && iter % 20 == 0) {
            solver.increaseBeta();
            std::cout << "Iter " << iter << ": Obj=" << obj << " (Beta updated)" << std::endl;
            solver.printTopology();
        }
        else if (iter % 10 == 0) {
            std::cout << "Iter " << iter << ": Obj=" << obj << std::endl;
            //solver.printTopology();
        }

        double max_change = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            double diff = std::abs(x[i] - x_old[i]);
            if (diff > max_change) max_change = diff;
        }
        if (iter % 10 == 0) {
            std::cout << "Iter " << iter << ": Max Change = " << max_change << std::endl;
        }

        // 如果变化极其微小，且已经迭代了一定次数（防止初期假收敛），则停止
        if (max_change < change_tolerance && iter > 100) {
            std::cout << ">>> CONVERGED at iter " << iter << "! Max change < " << change_tolerance << std::endl;
            break; // 退出循环
        }
        // --- 新增：Step G - 强制体积修正 (Hard Limit) ---
        //double current_vol = 0.0;
        //for (double val : x) current_vol += val;
        //double current_frac = current_vol / (solver.nx * solver.ny);

        //// 如果体积超标，强制全局扣除一个值 (Global Shift)
        //if (current_frac > Vol_frac) {
        //    // 简单的比例控制算法：超得越多，扣得越狠
        //    // 这里的 0.5 是一个修正速率系数
        //    double shift = (current_frac - Vol_frac) * 0.5;

        //    std::cout << "Vol violation: " << current_frac << " -> shifting by " << shift << std::endl;

        //    for (double& val : x) {
        //        val -= shift;
        //        if (val < 0.0) val = 0.0; // 再次截断
        //        // 注意：这里不需要检查 >1.0，因为我们是往下减
        //    }
        //}
    }

    std::cout << "Optimization Finished." << std::endl;
    solver.printTopology();

    return 0;


  //  int max_iter = 200;
  //  for (int iter = 0; iter < max_iter; ++iter) {

  //      // --- Step A: 正向几何处理 (x -> rho) ---
  //      solver.applyFilter();      // 密度滤波
  //      solver.applyProjection();  // Heaviside 投影

  //      // --- Step B: 物理分析 (FEM) 与 灵敏度计算 ---
  //      // 在真实软件中，这里会调用 FEM 解算器计算 Ku=f
  //      // 这里我们模拟一个“伪物理场”：
  //      // 目标：连接左侧墙壁(x=0)和右下角(x=nx, y=ny)
  //      double J_diff = solver.computeObjectiveAndGradient();
		//std::cout << "J_diff: " << J_diff << std::endl;
  //      // 我们手动修改 sensitivity 来模拟物理需求：
  //      // 假设右下角受力，左侧固定，材料应该分布在这两点之间
  //      // 灵敏度越高，表示该位置越需要材料
  //      //std::vector<double>& sens = const_cast<std::vector<double>&>(solver.getSensitivities());
  //      //for (int i = 0; i < ny; ++i) {
  //      //    for (int j = 0; j < nx; ++j) {
  //      //        // 简单的启发式引导：距离对角线越近，灵敏度越高
  //      //        double dist_diag = std::abs((double)i / ny - (double)j / nx);
  //      //        sens[i * nx + j] = 1.0 / (dist_diag + 0.1);

  //      //        // 强制左侧墙壁高灵敏度 (模拟固定端)
  //      //        if (j < 2) sens[i * nx + j] = 10.0;
  //      //    }
  //      //}

  //      // --- Step C: 应用悬垂约束 ---
  //      // 随着迭代进行，逐渐增加悬垂约束的权重
  //      double overhang_weight = 50;// (iter > 10) ? 0.5 : 0.0;
  //      solver.applyOverhangConstraint(overhang_weight);

  //      // --- Step D: 反向链式法则 (dJ/drho -> dJ/dx) ---
  //      solver.applyChainRule();

  //      // --- Step E: 优化器更新 x ---
  //      // 获取处理完链式法则后的最终灵敏度 dJ/dx
  //      std::vector<double>& x = solver.getDesignVariables();
  //      const std::vector<double>& final_sens = solver.getSensitivities();

  //      // 执行 OC 更新
  //      optimalityCriteriaUpdate(x, final_sens, vol_frac);
  //      // --- Step F: 延拓法 (Continuation Scheme) ---
  //      // 每 10 步增加一次 beta，使边界变清晰
  //      if (iter % 5 == 0) {
  //          solver.increaseBeta();
  //          std::cout << "Iter " << iter << ": Updating Beta." << std::endl;
  //          solver.printTopology(); // 打印当前拓扑
  //      }
  //  }

  //  std::cout << "Optimization Finished." << std::endl;
  //  solver.printTopology();

  //  return 0;
}

