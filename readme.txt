TaihuStone/
├── limitstl/             # 免支撑优化py代码
├── model/                # 存放输入的 3D 模型原始文件
    ├── ori/              # 存放缩放后的模型，用来比较
├── para/                 # 存放配置文件或参数设置 (Parameters)
├── result/               # 运行结果输出目录（如优化后的数据）
├── src/                  # 源代码目录
    ├── FastNoiseLite.h            # 轻量级噪声生成库（常用于生成随机多孔结构或扰动）
    ├── fd_check_diagnostics.json  # 诊断/调试用的 JSON 配置文件，优化过程生成
    ├── globalPara.cpp / .h        # 全局参数管理，控制模拟与优化的物理/数值参数
    ├── gpu_topology_optimizer*.py # 基于 GPU 的拓扑优化核心脚本（包含备份与测试版本）
    ├── main.cpp                   # C++ 程序入口，负责调度算法流水线
    ├── modelGen.cpp / .h          # 模型生成模块，处理几何构建逻辑
    ├── MorseComplex.cpp / .h      # 莫尔斯复形（Morse Complex）算法实现，用于拓扑分析
    ├── pathQuery.cpp              # 路径查询或渗透性分析相关逻辑
    ├── resultComp.cpp / .h        # 莫尔斯复形（Morse Complex）算法实现——宽松版
    ├── sdf_out_tube.stl           # 中间生成的有向距离场 (SDF) 管道模型文件
    ├── selfSupVoxel.cpp / .h      # 自支撑（Self-supporting）体素化处理，针对 3D 打印优化
    ├── Tool.cpp / .h              # 通用工具函数集（数学计算、文件 IO 等）
    ├── visualize.py               # 基于 Python 的结果可视化脚本
    ├── volume_render.py           # 体渲染脚本（用于查看 SDF 或体素场）
    └── vs_compare.cmd             # 视觉对比或版本差异对比的批处理脚本
├── TaihuStone/           # 项目核心子模块或资源文件夹
├── x64/                  # 64位编译生成的二进制文件及中间件
├── .gitignore            # Git 忽略文件配置
├── control.bat           # 通过para中的配置单次运行程序的脚本
├── control_multi.bat     # 批量运行/多任务控制脚本
├── opt_model.npy         # 优化后的模型数据 (NumPy 格式)
├── origin_model.npy      # 原始模型数据 (NumPy 格式)
├── pullcode.bat          # 代码拉取/同步脚本
├── pushcode.bat          # 代码推送/备份脚本
├── readme.txt            # 项目说明文档
├── result.stl            # 最终生成的 STL 网格模型
├── result_change.stl     # 修改或迭代后的 STL 模型
├── TaihuStone.sln        # Visual Studio 解决方案文件
├── TaihuStone.vcxproj    # Visual Studio 项目文件
├── TaihuStone.vcxproj.*  # VS 项目用户配置及过滤项
├── test.npy / test2.npy  # 测试用的数据文件
└── visualize.bat         # 结果对比的可视化脚本


General results: modelGen.cpp #159: sample_interior_points
Run modes:
1. direct run 
2. run control/control_multi.bat with para.txt

Regular results: modelGen.cpp #159: sample_regular

