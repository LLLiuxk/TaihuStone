@echo off
chcp 65001 >nul
title 拓扑优化程序运行中...

echo ========================================
echo   拓扑优化程序启动
echo ========================================

echo 步骤1: 切换到虚拟环境目录...
cd /d "D:\VSprojects\TaihuStone\limitstl\taihu_env\Scripts"

echo 步骤2: 激活虚拟环境...
call activate

echo 步骤3: 切换到项目目录...
cd /d "D:\VSprojects\TaihuStone\limitstl"

echo 步骤4: 运行拓扑优化程序...
echo ========================================

python gpu_topology_optimizer.py

pause
