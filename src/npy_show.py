import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# 1. 加载数据
# 假设您的 npy 是 3D 数组 (D, H, W)
data = np.load('sdf_out_tube.npy')
# 这里我们用刚才生成的模拟数据代替
slice_data = data[100, :, :] # 选取 z=32 的切面
# slice_data = slice_data_raw - 0.5
# 2. 准备画布
fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

# 3. 设置归一化：强制 0 对应色图中心
# SDF 数据通常正负范围不对称，TwoSlopeNorm 是关键
divnorm = colors.TwoSlopeNorm(vmin=slice_data.min(), vcenter=0., vmax=slice_data.max())

# 4. 绘制热力图 (Background)
# 使用 'RdBu' (Red-Blue) 色图，alpha 降低一点让等值线更明显
im = ax.imshow(slice_data, cmap='RdBu', norm=divnorm, origin='lower', alpha=0.9)

# 5. 绘制等值线 (Contours)
# (A) 辅助等值线：细、半透明
levels = np.linspace(slice_data.min(), slice_data.max(), 15)
ax.contour(slice_data, levels=levels, colors='k', linewidths=0.5, alpha=0.2)

# (B) 零等值线 (Surface)：粗、深色、高亮
# 这是 SDF 的核心，必须一眼就能看到
cs_zero = ax.contour(slice_data, levels=[0], colors='#222222', linewidths=2.5)

# 6. 修饰图表 (Academic Style)
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Signed Distance', rotation=270, labelpad=15)
cbar.outline.set_linewidth(0) # 去掉 colorbar 边框更现代

ax.set_title(r"SDF Visualization $\phi(x) = 0$", fontsize=14, pad=10)
ax.axis('off') # 如果不需要坐标轴，直接关掉更干净

plt.show()