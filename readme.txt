TaihuStone/ 项目目录与说明文档
=========================================================

一、 项目目录结构
---------------------------------------------------------
├── limitstl/             # 免支撑优化相关的 Python 脚本目录
│   ├── filter_stl.py         # 提取 STL 模型中最大连通分量的处理脚本
│   ├── gpu_topology_optimizer.py # GPU 拓扑优化 Python 脚本
│   ├── requirements.txt      # Python 依赖项列表
│   ├── run_optimize.bat      # 运行优化任务的批处理脚本
│   └── taihu_env/            # 本地 Python 虚拟环境目录（已在 .gitignore 中忽略）
├── model/                # 存放输入的原始 3D STL 模型文件
├── para/                 # 存放配置文件和优化参数设置 (Parameters)
│   ├── comp/                 # 各种对比实验的参数配置（例如 namaqualand 等）
│   └── comp2/ ~ comp5/       # 其他不同批次的参数配置文件目录
├── result/               # 运行结果输出目录
│   ├── before sfo/           # 优化前的数据/模型
│   ├── after sfo/            # 优化后的数据/模型
│   ├── fabricated/           # 用于制造/打印的输出模型
│   ├── msc_result/           # 莫尔斯复形分析结果
│   └── para/                 # 运行结果对应的参数备份
├── src/                  # C++ 与 Python 核心源代码目录
│   ├── FastNoiseLite.h       # 轻量级噪声生成库（常用于生成随机多孔结构或扰动）
│   ├── MorseComplex.cpp / .h # 莫尔斯复形（Morse Complex）算法实现，用于拓扑分析
│   ├── Tool.cpp / .h         # 通用工具函数集（数学计算、矩阵变换、文件 IO 等）
│   ├── TpmsPipeline.cpp / .h # 三周期最小表面 (TPMS) 隐式场生成、骨架提取与渗透率分析管道
│   ├── fd_check_diagnostics.json # 优化过程中生成的诊断与调试配置文件
│   ├── globalPara.cpp / .h   # 全局参数管理，控制模拟与优化的物理/数值参数
│   ├── gpu_topology_optimizer.py # 基于 GPU 的拓扑优化核心算法脚本
│   ├── main.cpp              # C++ 程序入口，负责调度算法流水线与演示模式
│   ├── modelGen.cpp / .h     # 模型生成模块，处理多孔结构构建及几何体素化逻辑
│   ├── pathQuery.cpp         # 路径查询或渗透性分析相关逻辑
│   ├── resultComp.cpp / .h   # 莫尔斯复形（Morse Complex）算法的松弛版/对比版实现
│   ├── selfSupVoxel.cpp / .h # 自支撑（Self-supporting）体素化处理与分析，针对 3D 打印优化
│   ├── selfSupVoxel_clean.cpp # 清理与优化后的自支撑体素化算法实现
│   ├── visualize.py          # 基于 Python 的三维网格及体素结果可视化脚本
│   ├── volume_render.py      # 体渲染脚本（用于查看 SDF 场或体素密度分布）
│   └── vs_compare.cmd        # 视觉/结构对比的批处理辅助命令
├── TaihuStone/           # Visual Studio 生成的中间编译资产与项目工程文件
├── x64/                  # 64位 C++ 编译生成的二进制可执行文件及链接库
├── .gitignore            # Git 忽略文件配置（已配置过滤 VS 缓存及本地 Python 环境）
├── control.bat           # 通过读取 para 中的配置单次运行主程序的脚本
├── control_compare.bat   # 运行模型对比分析的控制批处理脚本
├── glb_to_stl.py         # 将 GLB 格式模型批量或单文件转换为 STL 格式的 Python 脚本
├── pullcode.bat          # 代码拉取与同步脚本
├── pushcode.bat          # 代码推送与同步备份脚本
├── readme.txt            # 本项目说明文档（当前文档）
├── TaihuStone.sln        # Visual Studio 解决方案文件
├── TaihuStone.vcxproj    # Visual Studio 项目文件
├── TaihuStone.vcxproj.filters # Visual Studio 项目过滤器文件
├── visualize.bat         # 结果对比的可视化运行脚本
├── simplify.bat          # STL 模型自动简化入口，在 model 目录下运行，批量将超过 5MB 的模型简化至 5MB 以下
└── simplify_stls.py      # STL 简化 Python 脚本，基于 trimesh 和二次误差度量（QEM）算法进行简化


二、 运行模式与主要功能
---------------------------------------------------------
C++ 程序在 main.cpp 中支持两种运行模式 (通过修改 main.cpp 中的 chosen_fun 进行切换):
1. **chosen_fun = 0 (多孔模型生成与自支撑优化)**:
   - 读取指定的参数文件（如 `namaqualand_128_60_default_r.txt`），加载 STL 模型。
   - 对模型进行多孔结构化设计，检测自支撑性，输出优化后的多孔 STL 结果。
   - 运行方式: 双击运行 `control.bat`（默认加载 namaqualand 的参数文件）。

2. **chosen_fun = 1 (TPMS 骨架图分析与管道预览)**:
   - 生成常规及受扰动的 TPMS (Gyroid) 隐式场。
   - 对其进行体素化和拓架化 (Thinning)，生成骨架线并构建拓扑图。
   - 计算图节点的视觉渗透率 (VP, Visual Permeability)，最后将表面与骨架导出为 STL 模型，保存在 `tpms_output/` 下。

3. **模型简化功能 (simplify.bat)**:
   - 目的: 在进行大文件传输或网络同步时，原模型可能过大。本工具自动检测当前目录下所有大小超过 5MB 的 `.stl` 文件。
     * 简化后的模型统一输出存放在当前目录的 `simp_models/` 文件夹下。
   - 运行方式: 直接在当前目录下双击运行 `simplify.bat`。


三、 运行前依赖安装说明
---------------------------------------------------------
项目中的 Python 脚本（可视化、网格简化、格式转换等）需要以下第三方库：
- 可视化与处理: `trimesh`, `scipy`, `numpy`, `matplotlib`
- 二维简化内核: `fast-simplification` (由 `simplify_stls.py` 使用)

可以使用以下命令快速安装上述库：
  pip install trimesh numpy scipy matplotlib fast-simplification
