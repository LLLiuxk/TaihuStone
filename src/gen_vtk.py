import pyvista as pv
import numpy as np

# 1. 读取STL文件
mesh = pv.read("D:/VSprojects/TaihuStone/model/namaqualand.stl")

# 2. 创建体素网格（体素化）
# 设置体素大小（值越小分辨率越高）
voxel_size = 1.0  # 根据模型尺寸调整

# 创建体素网格
# 2. 创建体素网格（体素化）
# 使用新的API调用方式
voxel_size = 1.0  # 根据模型尺寸调整

# 修复点1：使用新的voxelize API
voxels = mesh.voxelize(density=voxel_size, check_surface=False)

# 3. 添加标量场用于体渲染
# 方法A：二值场（内部=1，外部=0）
voxels['scalar_field'] = np.ones(voxels.n_cells)

# 修复点2：正确计算隐式距离场
# 首先需要创建包含所有点的点云网格
if hasattr(voxels, 'cell_centers'):
    # 获取体素中心点
    centers = voxels.cell_centers().points
else:
    # 或者使用体素的点
    centers = voxels.points

# 计算体素中心到表面的距离
signed_dist = mesh.compute_implicit_distance(centers)
voxels['distance_field'] = signed_dist['implicit_distance']

# 4. 保存为VTI格式
voxels.save("D:/VSprojects/TaihuStone/model/output.vti")

print(f"体素网格信息：")
print(f"  体素数：{voxels.n_cells}")
print(f"  网格尺寸：{voxels.dimensions}")

