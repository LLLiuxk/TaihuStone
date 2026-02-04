"""
GPU加速的拓扑优化实现
使用CuPy和PyTorch进行GPU加速计算
根据overhang_constraints.tex文档实现
"""
import argparse  # <--- 1. 引入参数解析库
import sys
import numpy as np
import os
import torch
import torch.nn.functional as F
from scipy.sparse import csc_matrix
from skimage import measure
import trimesh
import time

# 检查GPU可用性
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"使用设备: {device}")

class GPUMMAOptimizer:
    """GPU加速的标准化 MMA 子问题求解（p/q 渐近线 + 原始闭式 + 对偶内循环）。

    接口保持不变：update(x, f0val, df0dx, fval, dfdx, xmin, xmax)
    - 渐近线 L/U 由 move_limit 维护
    - 子问题采用 Svanberg 形式的可分离保守近似：
        f̃0 = Σ (p0/(x-L) + q0/(U-x)) + const
        g̃i = Σ (Pi/(x-L) + Qi/(U-x)) + ri ≤ 0
      其中 p0,q0,Pi,Qi 由梯度和渐近线构造以匹配一阶导数。
    - 原始变量闭式：x(λ) = (√p̄·U + √q̄·L)/(√p̄ + √q̄)
    - 对偶变量 λ ≥ 0 通过投影牛顿/梯度法求解 g̃(x(λ)) ≤ 0。
    """

    def __init__(self, n_vars, move_limit=0.3, device='cuda'):
        self.n = n_vars
        self.move_limit = move_limit
        self.device = device

    def update(self, x, f0val, df0dx, fval, dfdx, xmin, xmax):
        """简化的MMA风格更新：
        - 以目标梯度为主方向
        - 对正的约束违反按幅度加权叠加梯度
        - 应用 move_limit 与 [xmin, xmax] 投影
        """
        # 组合步长方向：目标 + 违反的约束
        step = -df0dx.clone()
        # 违反量（仅对 >0 的约束施加惩罚）
        viol = torch.clamp(fval, min=0.0)
        if torch.any(viol > 0):
            weights = viol / (viol.sum() + 1e-8)
            # 减去约束梯度（因为我们要减小约束值）
            constraint_grad = dfdx.t().mm(weights.view(-1, 1)).view_as(step)
            # 增加惩罚力度：如果违反严重，加大约束梯度的权重
            # 显著增加惩罚系数以强制满足约束
            # 再次回调惩罚力度，接近"很好"的状态
            penalty_factor = 2000.0 + 20000.0 * torch.max(viol).item()
            #penalty_factor = 220.0 + 2200.0 * torch.max(viol).item()
            step -= penalty_factor * constraint_grad
        # 限制每步移动
        step = torch.clamp(step, -self.move_limit, self.move_limit)
        x_new = x + step
        # 投影到边界
        x_new = torch.max(torch.min(x_new, xmax), xmin)
        return x_new

class GPUTopologyOptimizer:
    """GPU加速的拓扑优化器（支持大模型）"""
    
    #def __init__(self, print_direction=[0, 0, 1], filter_radius=3.5, device='cuda', use_chunking=False):
    def __init__(self, print_direction=[0, 0, 1], filter_radius=3, device='cuda', use_chunking=False):
        self.device = device
        self.use_chunking = use_chunking
        self.chunk_size = 64 if use_chunking else None
        self.print_direction = torch.tensor(print_direction, dtype=torch.float32, device=device)
        self.print_direction = self.print_direction / torch.norm(self.print_direction)
        self.filter_radius = filter_radius
        
        # 优化参数 - 进一步调整以避免棋盘格模式
        self.beta_init = 1.0      # 初始值
        self.beta_max = 15.0      # 降低最大值，避免过度锐化, pre:12
        self.alpha_init = 1.0     
        self.alpha_max = 20.0     
        self.theta_crit = np.pi/4 
        self.t_support = 0.5      
        self.eps = 1e-6
        
        # MMA优化器
        self.mma = None
        
        if use_chunking:
            print(f"启用分块处理，块大小: {self.chunk_size}")

        # ROI（最小实体包围盒）信息
        self.roi_slices = None  # (slice_x, slice_y, slice_z)
        self.original_shape = None
        
    def create_base_protection_mask(self, target, safe_height=None):
        """创建底部保护掩码，用于屏蔽底面的悬挑约束"""
        if safe_height is None:
            safe_height = max(5, int(self.filter_radius * 2.5))
            
        # 确定主打印方向
        abs_dir = torch.abs(self.print_direction)
        main_axis = torch.argmax(abs_dir).item()
        
        # 置换使打印轴为 Z (dim=2)
        if main_axis == 0: # X
            target_perm = target.permute(1, 2, 0)
        elif main_axis == 1: # Y
            target_perm = target.permute(0, 2, 1)
        else:
            target_perm = target

        # target: (nx, ny, nz) -> (u, v, w)
        # 确保是二值掩码
        is_solid = (target_perm > 0.5).float()
        
        # 1. 找到每列是否有实体
        col_has_solid = (torch.sum(is_solid, dim=2) > 0.5) # (u, v)
        
        # 2. 找到每列第一个实体的位置 (W轴)
        # argmax 返回第一个最大值(1)的索引
        z_bottom = torch.argmax(is_solid, dim=2) # (u, v)
        
        # 3. 生成 mask
        nx, ny, nz = target_perm.shape
        z_indices = torch.arange(nz, device=self.device).view(1, 1, -1) # (1, 1, nz)
        z_bottom_expanded = z_bottom.unsqueeze(-1) # (nx, ny, 1)
        
        # 屏蔽条件：z < z_bottom + safe_height
        # 我们希望 mask=0 表示屏蔽，mask=1 表示保留
        # 所以当 z >= z_bottom + safe_height 时 mask=1
        # 另外，如果该列没有实体，mask=1 (不屏蔽，虽然也没东西)
        protection_mask = (z_indices >= (z_bottom_expanded + safe_height))
        
        # 对于没有实体的列，保持为 1 (True)
        protection_mask = protection_mask | (~col_has_solid.unsqueeze(-1))
        
        # 恢复维度
        if main_axis == 0:
            protection_mask = protection_mask.permute(2, 0, 1)
        elif main_axis == 1:
            protection_mask = protection_mask.permute(0, 2, 1)
        
        return protection_mask.float()

    def process_in_chunks(self, x, operation, *args, **kwargs):
        """分块处理大张量"""
        if not self.use_chunking:
            return operation(x, *args, **kwargs)
        
        nx, ny, nz = x.shape
        chunk_size = self.chunk_size
        results = []
        
        for i in range(0, nx, chunk_size):
            for j in range(0, ny, chunk_size):
                for k in range(0, nz, chunk_size):
                    # 提取块
                    chunk = x[i:min(i+chunk_size, nx), 
                             j:min(j+chunk_size, ny), 
                             k:min(k+chunk_size, nz)]
                    
                    # 处理块
                    chunk_result = operation(chunk, *args, **kwargs)
                    results.append((i, j, k, chunk_result))
        
        # 重组结果
        result_shape = results[0][3].shape if len(results) > 0 else x.shape
        if len(result_shape) == 3:
            full_result = torch.zeros(x.shape, device=self.device)
        else:
            full_result = torch.zeros((*x.shape, result_shape[-1]), device=self.device)
        
        for i, j, k, chunk_result in results:
            full_result[i:i+chunk_result.shape[0], 
                       j:j+chunk_result.shape[1], 
                       k:k+chunk_result.shape[2]] = chunk_result
        
        return full_result
        
    def setup_filter_gpu(self, shape):
        """设置GPU滤波器 - 标准高斯滤波"""
        self.shape = shape
        nx, ny, nz = shape
        
        # 创建标准高斯滤波器
        # 确保核大小为奇数
        kernel_size = max(int(4 * self.filter_radius + 1), 7)
        if kernel_size % 2 == 0:
            kernel_size += 1
        
        center = kernel_size // 2
        
        # 标准高斯核
        # sigma 取 filter_radius / 2.0 使滤波更平滑，抑制高频噪声
        sigma = self.filter_radius / 2.0
        kernel = torch.zeros((kernel_size, kernel_size, kernel_size), device=self.device)
        
        for i in range(kernel_size):
            for j in range(kernel_size):
                for k in range(kernel_size):
                    dist_sq = (i-center)**2 + (j-center)**2 + (k-center)**2
                    # 仅在半径范围内计算
                    if dist_sq <= (self.filter_radius * 2.0)**2:
                        kernel[i, j, k] = np.exp(-dist_sq / (2 * sigma**2))
        
        # 归一化
        kernel = kernel / torch.sum(kernel)
        
        self.filter_kernel = kernel.unsqueeze(0).unsqueeze(0)
        
    def filter_density_gpu(self, x, beta=None):
        """GPU密度滤波与投影"""
        # 1. 密度滤波 (卷积)
        x_batch = x.unsqueeze(0).unsqueeze(0)
        padding = self.filter_kernel.shape[-1] // 2
        # 使用 replicate 填充以避免边界处的密度下降
        x_padded = F.pad(x_batch, (padding, padding, padding, padding, padding, padding), mode='replicate')
        # pad mode='replicate' 意味着边界值向外复制，
        # 例如 [1, 1, 0] -> [1, 1, 1, 0, 0] (如果padding=1, 右边复制0? 不，右边复制最右边的值0，左边复制1)
        # 如果内部是 [1, ..., 1]，那么填充后全是1，卷积结果也是1。此处使用 valid 卷积 (padding=0)，因为已经手动填充了
        rho_filtered = F.conv3d(x_padded, self.filter_kernel, padding=0)
        rho_filtered = rho_filtered.squeeze(0).squeeze(0)
    
        
        # 2. Heaviside投影
        # 使用传入的 beta 或默认值
        if beta is None:
            p_beta = 1.0 # 默认值降低，初始更平滑
        else:
            p_beta = beta
            
        # 投影：将 [0,1] 映射到 [0,1]，但在 0.5 附近变陡
        # rho_projected = (tanh(beta * 0.5) + tanh(beta * (rho - 0.5))) / (2 * tanh(beta * 0.5))
        
        # 计算 tanh(beta * 0.5)，处理 beta 为 float 或 Tensor 的情况
        if isinstance(p_beta, (float, int)):
            tanh_beta_05 = np.tanh(p_beta * 0.5)
        else:
            tanh_beta_05 = torch.tanh(p_beta * 0.5)
            
        num = tanh_beta_05 + torch.tanh(p_beta * (rho_filtered - 0.5))
        den = 2 * tanh_beta_05
        rho_projected = num / den
        
        return rho_projected
    
    def compute_gradient_gpu(self, rho_tilde):
        """GPU计算密度梯度"""
        # 使用torch.gradient计算梯度
        grad_x = torch.gradient(rho_tilde, dim=0)[0]
        grad_y = torch.gradient(rho_tilde, dim=1)[0] 
        grad_z = torch.gradient(rho_tilde, dim=2)[0]
        
        return torch.stack([grad_x, grad_y, grad_z], dim=-1)

    def _compute_bbox_from_mask(self, mask_np, thr=0.5):
        """从二值或概率掩码计算最小包围盒（含至少一个体素）。返回 ((x0,x1),(y0,y1),(z0,z1))。"""
        idx = np.argwhere(mask_np > thr)
        if idx.size == 0:
            # 无实体，退化为全域
            sx, sy, sz = mask_np.shape
            return (0, sx), (0, sy), (0, sz)
        (x0, y0, z0) = idx.min(axis=0)
        (x1, y1, z1) = idx.max(axis=0)
        # 闭区间转半开区间
        return (int(x0), int(x1) + 1), (int(y0), int(y1) + 1), (int(z0), int(z1) + 1)

    def _expand_bbox(self, bbox, shape, margin):
        """对 bbox 以 margin 做各向同性扩展并裁剪到 shape 内。"""
        (sx, sy, sz) = shape
        (x0, x1), (y0, y1), (z0, z1) = bbox

        # 确定主打印方向，不扩展底板方向的 Margin，防止底部悬空
        # 假设打印方向为轴正向，且 AM filter 总是从 index 0 开始生长
        abs_dir = torch.abs(self.print_direction)
        main_axis = torch.argmax(abs_dir).item()
        
        # 如果主轴是X (0)，则 x0 (底部) 不扩展
        mx0 = 0 if main_axis == 0 else margin
        # 如果主轴是Y (1)，则 y0 (底部) 不扩展
        my0 = 0 if main_axis == 1 else margin
        # 如果主轴是Z (2)，则 z0 (底部) 不扩展
        mz0 = 0 if main_axis == 2 else margin

        x0 = max(0, x0 - mx0); x1 = min(sx, x1 + margin)
        y0 = max(0, y0 - my0); y1 = min(sy, y1 + margin)
        z0 = max(0, z0 - mz0); z1 = min(sz, z1 + margin)
        
        return (x0, x1), (y0, y1), (z0, z1)

    def _apply_roi_crop(self, x_np, target_np, margin):
        """根据目标掩码计算最小实体包围盒并扩展 margin，返回裁剪后的 (x_crop, target_crop) 与保存切片。"""
        bbox = self._compute_bbox_from_mask(target_np, thr=0.5)
        bbox_exp = self._expand_bbox(bbox, target_np.shape, margin)
        (x0, x1), (y0, y1), (z0, z1) = bbox_exp
        sx = slice(x0, x1); sy = slice(y0, y1); sz = slice(z0, z1)
        self.roi_slices = (sx, sy, sz)
        self.original_shape = target_np.shape
        return x_np[sx, sy, sz].copy(), target_np[sx, sy, sz].copy()
    
    def overhang_constraint_gpu(self, rho_tilde, gradient, beta, base_mask=None):
        """GPU悬挑角约束：修正版实现
        约束条件：(cos_angle - cos_critical) * grad_magnitude - small_positive ≤ 0
        """
        # 计算梯度幅值和方向
        grad_norm = torch.norm(gradient, dim=-1) + self.eps
        
        # 计算密度梯度方向与打印方向的余弦值
        # 注意：这里用正梯度方向，表示密度增加最快的方向
        cos_angle = torch.sum(self.print_direction * gradient, dim=-1) / grad_norm
        
        # 临界角的余弦值
        cos_critical = np.cos(self.theta_crit)

        # 修正的约束函数：(cos_angle - cos_critical) * grad_magnitude - threshold
        # 当这个值 > 0 时表示违反约束
        threshold = 0.0001  # 降低阈值，更严格
        constraint_value = (cos_angle - cos_critical) * grad_norm - threshold
        
        # 使用 ReLU 作为平滑的惩罚函数，替代硬截断
        # 这样当违反约束时，梯度可以回传到 cos_angle (即梯度方向)，允许优化器旋转法线
        # 而不仅仅是降低密度
        penalty_value = F.relu(constraint_value)
        
        # 用密度加权，只在有材料的地方施加约束
        # 使用平方项增加对大违反的惩罚力度 (Quadratic Penalty)
        # 这样对大的违反会有更大的梯度，迫使优化器优先解决严重区域
        #phi_ov = (penalty_value ** 2) * rho_tilde
        phi_ov = penalty_value * rho_tilde * 50.0  # 再次增强系数，确保在细节保留权重提高时仍能切除倒刺
        # --- 底部基座保护 ---
        # 使用传入的 base_mask 屏蔽底面约束
        if base_mask is not None:
            phi_ov = phi_ov * base_mask
        
        return phi_ov
    
    def hanging_constraint_gpu(self, rho_tilde, beta):
        """GPU悬挂特征约束：严格按照论文公式实现
        φ_i^(hg) = H_β(t ρ̃_i - S_i)
        其中 S_i = Σ_{k∈S(i)} w_{ik}^(s) ρ̃_k
        同时考虑“地面”支撑：在支撑方向的边界之外以常量 1 填充。
        """
        nx, ny, nz = rho_tilde.shape
        t = self.t_support  # 支撑阈值参数 (0.3~0.6)

        # 找到主轴（打印方向）与支撑方向符号（与打印方向相反）
        abs_dir = torch.abs(self.print_direction)
        main_axis = torch.argmax(abs_dir).item()
        sgn = float((-torch.sign(self.print_direction[main_axis])).item())

        # 创建支撑区域卷积核 S(i)
        kernel_size = 5  # 扩大支撑区域
        support_kernel = torch.zeros((1, 1, kernel_size, kernel_size, kernel_size), device=self.device)
        center = kernel_size // 2

        # 设置支撑区域权重（沿打印方向的下方）
        for i in range(kernel_size):
            for j in range(kernel_size):
                for k in range(kernel_size):
                    pos = [i - center, j - center, k - center]
                    if main_axis == 0 and pos[0] * (-torch.sign(self.print_direction[0])) > 0:
                        dist = np.sqrt(pos[1] ** 2 + pos[2] ** 2)
                        support_kernel[0, 0, i, j, k] = max(0.0, 2 - dist)
                    elif main_axis == 1 and pos[1] * (-torch.sign(self.print_direction[1])) > 0:
                        dist = np.sqrt(pos[0] ** 2 + pos[2] ** 2)
                        support_kernel[0, 0, i, j, k] = max(0.0, 2 - dist)
                    elif main_axis == 2 and pos[2] * (-torch.sign(self.print_direction[2])) > 0:
                        dist = np.sqrt(pos[0] ** 2 + pos[1] ** 2)
                        support_kernel[0, 0, i, j, k] = max(0.0, 2 - dist)

        # 归一化权重
        if torch.sum(support_kernel) > 0:
            support_kernel = support_kernel / torch.sum(support_kernel)

        # 计算支撑量度 S_i（带“地面”支撑）：在支撑方向外部以 1 填充
        rho_batch = rho_tilde.unsqueeze(0).unsqueeze(0)  # [1,1,D,H,W]
        pad = center
        # F.pad 的顺序为 (W_left, W_right, H_left, H_right, D_left, D_right)
        xpad = F.pad(rho_batch, (pad, pad, pad, pad, pad, pad), mode='constant', value=0.0)
        # 仅在支撑侧填充 1 层“地面”厚度
        ground_thickness = 1
        if main_axis == 0:
            if sgn < 0:
                xpad[:, :, :ground_thickness, :, :] = 1.0
            else:
                xpad[:, :, -ground_thickness:, :, :] = 1.0
        elif main_axis == 1:
            if sgn < 0:
                xpad[:, :, :, :ground_thickness, :] = 1.0
            else:
                xpad[:, :, :, -ground_thickness:, :] = 1.0
        else:  # main_axis == 2
            if sgn < 0:
                xpad[:, :, :, :, :ground_thickness] = 1.0
            else:
                xpad[:, :, :, :, -ground_thickness:] = 1.0

        # 卷积（无额外 padding），输出大小与原张量一致
        S_i = F.conv3d(xpad, support_kernel, padding=0).squeeze(0).squeeze(0)

        # 论文公式：v_i = t ρ̃_i - S_i
        v_i = t * rho_tilde - S_i

        # 光滑Heaviside函数：φ_i^(hg) = H_β(v_i)
        phi_hg = self._heaviside_gpu(v_i, beta)
        return phi_hg
    
    def _heaviside_gpu(self, u, beta):
        """GPU光滑Heaviside函数
        
        理想的Heaviside函数：
        - u > 0 时，H(u) = 1
        - u ≤ 0 时，H(u) = 0
        
        使用光滑近似：H_β(u) = 1 / (1 + e^(-β*u))
        - 当β很大时，函数接近理想的阶跃函数
        - u > 0 时，H_β(u) → 1
        - u < 0 时，H_β(u) → 0
        - u = 0 时，H_β(u) = 0.5
        """
        return 1.0 / (1.0 + torch.exp(-beta * u))
    
    def _sigmoid_gpu(self, rho, beta):
        """GPU光滑sigmoid函数"""
        return rho
    
    def aggregate_constraint_gpu(self, phi):
        """GPU约束聚合：直接求和聚合
        G(x) = Σ_i φ_i ≤ 0
        
        根据图片中的公式，约束聚合应该是所有局部违反量的直接求和，
        而不是log-sum-exp形式。这更符合传统的约束处理方法。
        """
        # 直接对所有违反量求和
        G = torch.sum(phi)
        
        return G
    
    def apply_am_filter(self, rho):
        """应用增材制造(AM)滤波器模拟逐层打印过程
        
        原理：
        每一层的材料必须由下一层的材料支撑。
        rho_printed[z] = min(rho[z], dilate(rho_printed[z-1]))
        
        改进：使用 SoftMin 替代 HardMin，确保梯度能够同时流向 rho 和 support。
        这解决了初始化均匀时梯度无法向下传播导致支撑结构无法生长的问题。
        """
        # 确定主打印方向
        abs_dir = torch.abs(self.print_direction)
        main_axis = torch.argmax(abs_dir).item()
        
        # 如果不是 Z 轴 (2)，则置换维度使打印轴为最后一维
        # 原始: (nx, ny, nz)
        if main_axis == 0: # X轴打印
            # 变为 (ny, nz, nx)，此时 nx 是层
            rho_working = rho.permute(1, 2, 0) 
        elif main_axis == 1: # Y轴打印
            # 变为 (nx, nz, ny)，此时 ny 是层
            rho_working = rho.permute(0, 2, 1)
        else:
            rho_working = rho
            
        nx, ny, nz = rho_working.shape
        
        # 使用列表收集每一层
        layers = []
        
        # 第0层（底层）假设完全由基板支撑
        layers.append(rho_working[:, :, 0])
        
        # SoftMin 参数
        # smin(a, b) = (a + b - sqrt((a-b)^2 + eps)) / 2
        eps = 1e-6
        
        for k in range(1, nz):
            # 获取上一层的打印状态 (nx, ny) -> (1, 1, nx, ny)
            prev_layer = layers[-1].unsqueeze(0).unsqueeze(0)
            
            # 计算上一层能提供的支撑范围 (3x3 max pool)
            support = F.max_pool2d(prev_layer, kernel_size=3, stride=1, padding=1)
            support = support.squeeze(0).squeeze(0)
            
            # SoftMin: 平滑的最小值函数
            # 允许梯度在 a=b 附近同时流向两者
            a = rho_working[:, :, k]
            b = support
            current_layer = (a + b - torch.sqrt((a - b)**2 + eps)) / 2.0
            
            layers.append(current_layer)
            
        # 堆叠所有层
        rho_printed = torch.stack(layers, dim=2)
        
        # 恢复维度
        if main_axis == 0:
            rho_printed = rho_printed.permute(2, 0, 1)
        elif main_axis == 1:
            rho_printed = rho_printed.permute(0, 2, 1)
            
        return rho_printed

    def compute_gradients_gpu(self, x, rho_tilde, gradient, beta, alpha, target, base_mask=None):
        """GPU计算所有梯度 - 基于AM滤波器的新版"""
        # 启用梯度计算
        x.requires_grad_(True)
        
        # 1. 物理密度滤波 (消除棋盘格，控制最小尺寸)
        # 传入 beta 以控制投影锐度
        rho_phys = self.filter_density_gpu(x, beta=beta)
        
        # 2. AM 打印过程模拟 (处理自支撑约束)
        # rho_printed 代表了"实际能打印出来的结构"
        rho_printed = self.apply_am_filter(rho_phys)
        
        # 3. 目标函数优化
        # 使用"结构感知"的目标函数：
        # 1. 锐化目标值：使用sigmoid将连续的目标值在0.5附近拉开，
        #    使得 0.51 -> ~1.0, 0.49 -> ~0.0，正确反映其拓扑属性。
        # 2. 结合 Soft Dice Loss：最大化重叠区域，比L2更能捕捉几何形状。
        # 3. 保留加权MSE：但针对锐化后的目标场计算。
        
        sharpness = 100.0
        # 将目标场转换为实际上我们想要匹配的"结构场"
        target_structure = torch.sigmoid((target - 0.5) * sharpness)
        
        # A. Soft Dice Loss (最大化结构重叠，范围 [0, 1])
        # Dice = 2*|A∩B| / (|A| + |B|)
        intersection = torch.sum(rho_printed * target_structure)
        union = torch.sum(rho_printed) + torch.sum(target_structure) + 1e-6
        # 我们要最小化 Loss，所以用 1 - Dice
        # 放大 Dice 权重使其梯度在数值上可观
        #dice_weight = 100000.0 
        dice_weight = 80000.0 
        loss_dice = (1.0 - (2.0 * intersection / union)) * dice_weight
        
        # B. Weighted MSE (像素级强度匹配)
        # 权重矩阵：前景区域权重更高
        weights = 1.0 + 9.0 * target_structure
        diff = rho_printed - target_structure
        # 使用 sum 而不是 mean，保持与其他项量级一致
        loss_mse = torch.sum(weights * (diff**2))

        #f0 = loss_dice + loss_mse
          # C. Total Variation (TV) Regularization - 抑制棋盘格而不模糊边缘
        # TV = sum(|grad_x| + |grad_y| + |grad_z|)
        # 这种正则化偏好分块常数解(Piecewise Constant)，能保留锐利边缘，同时去除高频振荡
        rho_p = rho_printed
        dx = torch.abs(rho_p[1:,:,:] - rho_p[:-1,:,:])
        dy = torch.abs(rho_p[:,1:,:] - rho_p[:,:-1,:])
        dz = torch.abs(rho_p[:,:,1:] - rho_p[:,:,:-1])
        tv_loss = torch.sum(dx) + torch.sum(dy) + torch.sum(dz)
        
        # 权重设置：TV 不需要太大，只需足够抑制高频噪声
        # 降低 TV 权重以允许合法的表面纹理存在
        tv_weight = 0.1 
        loss_tv = tv_weight * tv_loss

        f0 = loss_dice + loss_mse + loss_tv
        
        # 计算显式悬挑约束值用于监控（不参与梯度计算）
        # 真正的约束是通过 AM 滤波器隐式实现的
        # 这里我们计算 rho_phys 和 rho_printed 的差异作为"不可打印体积"的度量
        
        # 差异越大，说明 rho_phys 中有越多部分被 AM 滤波器"削减"掉了
        # 即这些部分是悬空的
        unprintable_volume = torch.sum(torch.abs(rho_phys - rho_printed))
        G_ov = unprintable_volume.item()
        
        # --- 策略变更 ---
        # 仅仅使用惩罚函数(Soft Constraint)会导致最终结果仍有残留。
        # 我们现在将 G_ov 作为一个显式的硬约束传递给优化器。
        # 为了优化器稳定性，我们在这里不把 punishment 加到 f0 中，而是分别计算梯度。
        
        # 1. 目标函数梯度 (仅 Loss + MSE)
        f0.backward(retain_graph=True)
        df0dx = x.grad.clone()
        x.grad.zero_()
        
        # 2. 悬挑约束梯度 (unprintable_volume)
        # 这是一个显式约束： unprintable_volume <= tolerance
        unprintable_volume.backward()
        dg_ov = x.grad.clone()
        x.grad.zero_() # 清理梯度
        
        G_hg = 0.0
        dg_hg = torch.zeros_like(df0dx)
        
        x.requires_grad_(False)
        
        return f0.item(), df0dx, G_ov, dg_ov, G_hg, dg_hg
    
    def optimize(self, initial_x, target, output_dir, max_iter=50, use_roi=False):
        """GPU主优化循环（支持大模型）

        说明:
        - 每次迭代导出中间 STL/NPY 到 iter_outputs/
        - 使用连续密度（level=0.5）进行等值面提取
        - 仅在目标的最小实体包围盒（加 margin）内优化，避免边界“板子”
        """
        print("开始GPU加速拓扑优化")
        print(f"设备: {self.device}")
        print(f"模型大小: {initial_x.shape}")
        print("- 并行计算")
        print("- GPU内存优化")
        print("- 自动梯度计算")

        if use_roi:
            # ROI 裁剪（基于目标掩码的最小实体包围盒，带 margin）
            margin = 2
            x_np_crop, target_np_crop = self._apply_roi_crop(initial_x, target, margin)
            print(f"ROI 裁剪: 原始 {initial_x.shape} → ROI {x_np_crop.shape}，margin={margin}")
            x = torch.tensor(x_np_crop, dtype=torch.float32, device=self.device)
            target_gpu = torch.tensor(target_np_crop, dtype=torch.float32, device=self.device)
        else:
            # 使用完整体素场，保持坐标与baseline对齐
            print("使用完整体素场优化，保持坐标对齐")
            self.roi_slices = None
            self.original_shape = initial_x.shape
            x = torch.tensor(initial_x, dtype=torch.float32, device=self.device)
            target_gpu = torch.tensor(target, dtype=torch.float32, device=self.device)

        # 设置滤波器（针对 ROI 尺寸）
        self.setup_filter_gpu(x.shape)

        # 记录 ROI 原点偏移（用于导出时对齐原坐标系）
        if self.roi_slices is not None:
            sx, sy, sz = self.roi_slices
            roi_origin = (sx.start or 0, sy.start or 0, sz.start or 0)
        else:
            roi_origin = (0, 0, 0)

        # 目标做同核预滤波，避免对齐偏置
        # [已移除] 目标滤波，改为直接使用原始 target
        # with torch.no_grad():
        #    target_filtered = self.filter_density_gpu(target_gpu)
            
        # 保存原始目标形状用于对比
        print("保存原始目标形状(original_target.stl)...")
        try:
            # 直接通过 NPY 数据生成 STL，不进行滤波
            if use_roi:
                self._export_stl_from_density(
                    target_gpu.cpu().numpy(),
                    output_dir + "/original_target.stl",
                    level=0.5,
                    step_size=1,
                    origin=roi_origin,
                )
            else:
                self._export_stl_from_density(
                    target_gpu.cpu().numpy(),
                    output_dir + "original_target.stl",
                    level=0.5,
                    step_size=1,
                )
        except Exception as e:
            print(f"保存原始目标STL失败: {e}")
            
        # 创建底部保护掩码（基于原始二值目标）
        # 这将屏蔽每列最底部的实体区域，防止其被误判为悬挑
        base_mask = self.create_base_protection_mask(target_gpu)
        print("已创建底部保护掩码，防止底面变形")

        # 初始化MMA
        # 进一步降低移动限制，确保生长过程平稳
        self.mma = GPUMMAOptimizer(x.numel(), move_limit=0.02, device=self.device)

        # 连续化参数
        beta = self.beta_init
        alpha = self.alpha_init

        history = []

        # 导出设置
        save_intermediate = True
        save_every = 5
        out_dir = "D:\\VSprojects\\TaihuStone\\limitstl\\iter_outputs"
        os.makedirs(out_dir, exist_ok=True)

        # 内存管理
        if self.device == 'cuda':
            torch.cuda.empty_cache()

        for iteration in range(max_iter):
            start_time = time.time()
            print(f"\n--- GPU拓扑优化迭代 {iteration+1}/{max_iter} ---")
            print(f"连续化参数: β={beta:.1f}, α={alpha:.1f}")

            # GPU内存监控
            if self.device == 'cuda':
                allocated = torch.cuda.memory_allocated() / 1e9
                reserved = torch.cuda.memory_reserved() / 1e9
                print(f"GPU内存: {allocated:.2f}/{reserved:.2f} GB")

            # 计算目标函数与约束及其梯度
            try:
                # 不再使用 filtering 后的 target，而是使用原始 target_gpu (已是float)
                # target_filtered 变量不再需要，直接传 target_gpu
                f0, df0dx, G_ov, dg_ov, G_hg, dg_hg = self.compute_gradients_gpu(
                    x, None, None, beta, alpha, target_gpu, base_mask=base_mask
                )
            except RuntimeError as e:
                if "out of memory" in str(e).lower():
                    print("GPU内存不足，释放缓存后重试下一轮…")
                    if self.device == 'cuda':
                        torch.cuda.empty_cache()
                    continue
                raise

            print(f"Obj: {f0:.6f}  Cons_overhang: {G_ov:.6f}")
            #print(f"悬挑约束: {G_ov:.6f} (目标: ≤0)")
            #print(f"悬挂约束: 已禁用")

            # MMA 更新 - 无显式约束
            # 动态调整移动限制：前期大步长加速生长，后期小步长精细化
            # 初始 0.15，衰减至 0.05，保持一定的搜索能力
            current_move_limit = max(0.005, 0.15 * (0.99 ** iteration))
            self.mma.move_limit = current_move_limit
            
            # 显式约束传递给 MMA
            # G_ov 是 total unprintable volume。
            # 引入容差 (Tolerance)，允许少量不可打印体积 (例如边缘效应)
            # 容差值设为 400 体素 (约占总体积 0.2%)，避免过度牺牲目标函数
            # now :10
            constraint_tolerance = 10.0
            scale_factor = 100.0
            
            # 构造约束向量: (G_ov - tolerance) / scale <= 0
            constraint_val = (G_ov - constraint_tolerance) / scale_factor
            
            fval = torch.tensor([constraint_val], device=self.device)
            # dfdx 应该是是一个 (m, n) 的 2D 矩阵
            # dg_ov 目前可能是 多维 Tensor，需要展平为 1D
            # 注意：约束梯度的方向不变，还是指向减小 G_ov 的方向
            dg_ov_flat = dg_ov.flatten()
            dfdx = (dg_ov_flat / scale_factor).unsqueeze(0)
            
            x_flat = x.flatten()
            x_new = self.mma.update(
                x_flat, f0, df0dx.flatten(), fval, dfdx,
                torch.zeros_like(x_flat), torch.ones_like(x_flat)
            )

            # 自适应阻尼：如果发生剧烈震荡（体素数量跳变），强制回退并减小步长
            # 这里简单地平滑步长，通过动量项减少高频震荡
            if iteration > 0:
                momentum = 0.5  # 增加动量系数以增强稳定性
                # 将新计算的x与上一次的x进行平均，抑制跳变
                x_new = (1 - momentum) * x_new + momentum * x_flat

            # 计算最大密度变化 (L-infinity norm)
            # 这比目标函数变化更能反映收敛状态
            x_new_reshaped = x_new.view_as(x)
            max_density_change = torch.max(torch.abs(x_new_reshaped - x)).item()

            x = x_new_reshaped # 更新 x

            # [已移除] 增强的平滑处理以减少棋盘格模式
            # 这里的额外滤波被移除，仅依赖 filter_density_gpu 进行必要的棋盘格控制

            # 记录与导出
            rho_current = self.filter_density_gpu(x)
            vol_current = torch.mean(rho_current).item()
            iter_time = time.time() - start_time
            history.append({
                'iteration': iteration,
                'objective': f0,
                'overhang_constraint': G_ov,
                'hanging_constraint': G_hg,
                'volume_fraction': vol_current,
                'max_change': max_density_change,
                'time': iter_time
            })
            
            print(f"迭代时间: {iter_time:.2f}s | 最大密度变化: {max_density_change:.6f}")
            if save_intermediate and (iteration % save_every == 0):
                base = os.path.join(out_dir, f"iter_{iteration:03d}")
                # [已禁用] 只保存STL，不再保存NPY以节省空间
                # try:
                #     if use_roi:
                #         # 保存 ROI 密度
                #         np.save(base + "_roi.npy", x.detach().cpu().numpy())
                #         # 回填到原始尺寸，方便全局对齐比较
                #         full = np.zeros(self.original_shape, dtype=np.float32)
                #         sx, sy, sz = self.roi_slices
                #         full[sx, sy, sz] = rho_current.detach().cpu().numpy()
                #         np.save(base + "_full.npy", full)
                #     else:
                #         # 直接保存完整密度，确保坐标对齐
                #         np.save(base + ".npy", rho_current.detach().cpu().numpy())
                # except Exception as _e:
                #     print(f"保存中间NPY失败: {_e}")
                try:
                    if use_roi:
                        # 直接导出 ROI STL，并对顶点添加原点偏移，保证与全局坐标一致
                        self._export_stl_from_density(
                            rho_current.detach().cpu().numpy(),
                            base + ".stl",
                            level=0.5,
                            step_size=2,
                            origin=roi_origin,
                        )
                    else:
                        # 直接导出完整STL，坐标与baseline对齐
                        self._export_stl_from_density(
                            rho_current.detach().cpu().numpy(),
                            base + ".stl",
                            level=0.5,
                            step_size=2,
                        )
                except Exception as _e:
                    print(f"保存中间STL失败: {_e}")

            # 收敛与连续化调整
            # 使用与MMA更新中一致的容差 (约400体素)
            constraint_tol = 10.0 
            constraints_satisfied = (G_ov <= constraint_tol)
            
            # 1. 密度场收敛判定 (最重要)
            # 如果设计变量不再发生显著变化，说明已经找到了(局部)极值
            is_density_converged = (max_density_change < 0.001)
            
            # 2. 目标函数平稳判定
            is_objective_stable = False
            if iteration > 10:
                 # 检查过去5次迭代的目标函数波动
                 recent_objs = [h['objective'] for h in history[-5:]]
                 obj_change = (max(recent_objs) - min(recent_objs)) / (abs(recent_objs[-1]) + 1e-6)
                 # 0.1% 的波动认为稳定
                 if obj_change < 0.001: 
                     is_objective_stable = True

            if is_density_converged:
                if constraints_satisfied:
                    print(f"收敛: 密度变化极小 ({max_density_change:.6f}) 且 满足悬挑约束 ({G_ov:.1f} <= {constraint_tol})")
                    break
                else:
                    print(f"密度已收敛 ({max_density_change:.6f}) 但未满足约束 ({G_ov:.1f} > {constraint_tol})。")
                    # 如果已经非常收敛但约束不满足，可能是参数问题，与其空转不如停止或调整beta
                    # 这里如果是多次连续收敛但违反约束，则强制停止
                    if iteration > 50 and all(h['max_change'] < 0.001 for h in history[-10:]):
                         print("无法进一步满足约束，强制停止优化。")
                         break

            if not constraints_satisfied:
                # 减缓 beta 增长速度，每 10 次迭代更新一次，且增长倍率降低，给予更多收敛时间
                if iteration > 0 and iteration % 10 == 0:
                    old_beta, old_alpha = beta, alpha
                    beta = min(self.beta_max, beta * 1.2)
                    alpha = min(self.alpha_max, alpha * 1.2)
                    if beta > old_beta or alpha > old_alpha:
                        print(f"连续化参数更新: β={old_beta:.1f}→{beta:.1f}, α={old_alpha:.1f}→{alpha:.1f}")
            else:
                # 约束满足情况下，如果目标函数也非常稳定，也可以停止
                if is_objective_stable and iteration > 20: 
                    print("悬挑约束满足且目标函数稳定")
                    break

            if (iteration % 5 == 0) and (self.device == 'cuda'):
                torch.cuda.empty_cache()

        return x.detach().cpu().numpy(), history

    def _export_stl_from_density(self, density, stl_path, level=0.5, step_size=1, origin=None):
        """使用连续密度场导出 STL（Marching Cubes）。

        参数:
          density: np.ndarray, 连续密度 [0,1]
          stl_path: 输出 STL 路径
          level: 等值面阈值（默认 0.5）
          step_size: MC 步长，>1 可显著减少三角形数量，适合中间导出
          origin: (x0,y0,z0) 顶点平移（用于 ROI 偏移对齐全局坐标）
        """
        field = np.asarray(density, dtype=np.float32)
        # 快速检查是否有足够体素超过阈值
        solid_voxels = int(np.count_nonzero(field >= level))
        #print(f"[导出] level={level}，体素>=level数量: {solid_voxels}")
        if solid_voxels < 50:
            print("Too less voxels, skip generating STL")
            return
        try:
            # 自动闭合网格：在四周填充一圈 0，强迫边界处生成封闭面
            # padding 1 pixel with 0
            pad_width = 1
            field_padded = np.pad(field, pad_width=pad_width, mode='constant', constant_values=0)
            
            verts, faces, normals, values = measure.marching_cubes(
                field_padded, level=level, spacing=(1.0, 1.0, 1.0), step_size=step_size, allow_degenerate=False
            )
            
            # 修正坐标：由于 padding 导致原点偏移了 (1,1,1)，需要减回去
            verts = verts - pad_width
            
            
            # ==========================================================
            # 关键修正 A：轴顺序调整 (解决角度问题)
            # ==========================================================
            # 目前 verts 是 (z, y, x) 顺序 (对应 numpy 的 axis 0, 1, 2)
            # STL 需要 (x, y, z)
            # 所以我们需要交换第 0 列和第 2 列
        
            #verts = verts[:, [2, 1, 0]]
            #res = density.shape[0]
            #voxel_size = 1.0 / (res - 1)

            # 2. 处理缩放逻辑
            #if voxel_size is not None:
                # 【修改点1】修正缩进，使其包含在 if voxel_size is not None 内部
                # 【修改点2】使用 np.isscalar 增强兼容性（防止 numpy float 类型报错）
                #if np.isscalar(voxel_size) or isinstance(voxel_size, (int, float)):
                 #   scale = np.array([voxel_size, voxel_size, voxel_size])
                #else:
                 #   scale = np.array(voxel_size) # 假设输入是 (sx, sy, sz)
                
                # 【修改点3】千万不要忘了应用缩放！
                #verts = verts * scale

            # 核心修正：如果C++里是居中的，那么体素网格的左下角起点应该是 -0.5
           # start_point = (-0.5, -0.5, -0.5)
            #origin = start_point

            if origin is not None:
                ox, oy, oz = origin
                verts = verts + np.array([ox, oy, oz], dtype=verts.dtype)
            if len(verts) == 0 or len(faces) == 0:
                print("Illegal mesh, skip!")
                return
            mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
            # 修复法向问题：翻转法向以指向外部
            # mesh.invert()
            mesh.export(stl_path)
            print(f"Intermediate STL: {stl_path}，V={len(verts)}, F={len(faces)}")
        except Exception as e:
            print(f"Write STL failed!: {e}")
    
    def finite_difference_check(self, x_numpy, target_numpy, n_checks=50, h=1e-4, beta=None, alpha=None):
        """对 G_ov 和 G_hg 做有限差分校验（中心差分），比较解析梯度（自动求导）和数值梯度。

        参数:
          x_numpy: numpy 数组，当前设计变量（密度场）
          target_numpy: numpy 数组，目标几何掩码
          n_checks: 随机检查的变量数量
          h: 差分步长
          beta, alpha: 连续化参数，若为 None 则使用初始化值
        返回: 字典，包含每个约束的误差统计
        """
        if beta is None:
            beta = self.beta_init
        if alpha is None:
            alpha = self.alpha_init

        device = self.device
        x_t = torch.tensor(x_numpy, dtype=torch.float32, device=device)
        target_t = torch.tensor(target_numpy, dtype=torch.float32, device=device)

        # 确保滤波器已创建并与输入尺寸匹配
        try:
            shape_matches = hasattr(self, 'shape') and tuple(self.shape) == tuple(x_t.shape)
        except Exception:
            shape_matches = False
        if (not hasattr(self, 'filter_kernel')) or (not shape_matches):
            # setup_filter_gpu 接受 python tuple 形状
            self.setup_filter_gpu(tuple(x_t.shape))

        # 计算解析梯度（使用GPU自动求导实现的函数）
        _, _, G_ov_val, dg_ov, G_hg_val, dg_hg = self.compute_gradients_gpu(
            x_t.clone(), None, None, beta, alpha, target_t
        )

        dg_ov_flat = dg_ov.flatten().cpu().numpy()
        dg_hg_flat = dg_hg.flatten().cpu().numpy()

        n_vars = x_t.numel()
        torch.manual_seed(0)
        idxs = torch.randperm(n_vars, device=device)[:min(n_checks, n_vars)].cpu().numpy()

        abs_err_ov = []
        rel_err_ov = []
        abs_err_hg = []
        rel_err_hg = []
        per_index_info = []

        # 中心差分评估每个选中变量
        for idx in idxs:
            # 构造 x_plus 和 x_minus（自适应局部步长）
            x_flat = x_t.clone().view(-1)
            current_val = float(x_flat[idx].item())
            h_local = float(h) * max(1.0, abs(current_val))

            x_plus = x_flat.clone()
            x_minus = x_flat.clone()
            x_plus[idx] = x_plus[idx] + h_local
            x_minus[idx] = x_minus[idx] - h_local
            x_plus = x_plus.view(x_t.shape)
            x_minus = x_minus.view(x_t.shape)

            # 使用 no_grad 前向计算 G_ov 和 G_hg
            with torch.no_grad():
                rho_p = self.filter_density_gpu(x_plus)
                grad_p = self.compute_gradient_gpu(rho_p)
                phi_ov_p = self.overhang_constraint_gpu(rho_p, grad_p, beta)
                G_ov_p = self.aggregate_constraint_gpu(phi_ov_p)
                # 暂时禁用悬挂约束
                # phi_hg_p = self.hanging_constraint_gpu(rho_p, beta)
                # G_hg_p = self.aggregate_constraint_gpu(phi_hg_p)
                G_hg_p = torch.tensor(0.0, device=self.device)

                rho_m = self.filter_density_gpu(x_minus)
                grad_m = self.compute_gradient_gpu(rho_m)
                phi_ov_m = self.overhang_constraint_gpu(rho_m, grad_m, beta)
                G_ov_m = self.aggregate_constraint_gpu(phi_ov_m)
                # 暂时禁用悬挂约束
                # phi_hg_m = self.hanging_constraint_gpu(rho_m, beta)
                # G_hg_m = self.aggregate_constraint_gpu(phi_hg_m)
                G_hg_m = torch.tensor(0.0, device=self.device)

            # 数值导数（使用自适应步长）
            num_dG_ov = (G_ov_p.item() - G_ov_m.item()) / (2.0 * h_local)
            num_dG_hg = (G_hg_p.item() - G_hg_m.item()) / (2.0 * h_local)

            ana_dG_ov = float(dg_ov_flat[idx])
            ana_dG_hg = float(dg_hg_flat[idx])

            abs_ov = abs(num_dG_ov - ana_dG_ov)
            abs_hg = abs(num_dG_hg - ana_dG_hg)

            # 更稳定的相对误差分母：使用 num 或 ana 的最大绝对值，且不低于一个阈值
            denom_ov = max(abs(num_dG_ov), abs(ana_dG_ov), 1e-8)
            denom_hg = max(abs(num_dG_hg), abs(ana_dG_hg), 1e-8)

            rel_ov = abs_ov / denom_ov
            rel_hg = abs_hg / denom_hg

            abs_err_ov.append(abs_ov)
            rel_err_ov.append(rel_ov)
            abs_err_hg.append(abs_hg)
            rel_err_hg.append(rel_hg)

            coord_tuple = np.unravel_index(int(idx), tuple(x_t.shape))
            coord_list = [int(ci) for ci in coord_tuple]
            per_index_info.append({
                'idx': int(idx),
                'coord': coord_list,
                'num_dG_ov': float(num_dG_ov),
                'ana_dG_ov': float(ana_dG_ov),
                'abs_err_ov': float(abs_ov),
                'num_dG_hg': float(num_dG_hg),
                'ana_dG_hg': float(ana_dG_hg),
                'abs_err_hg': float(abs_hg),
                'h_local': float(h_local)
            })

        import statistics
        result = {
            'n_checks': len(idxs),
            'G_ov': {
                'analytic_max_abs': float(max(map(abs, dg_ov_flat[idxs]))),
                'num_mean_abs_error': float(statistics.mean(abs_err_ov)) if len(abs_err_ov) > 0 else 0.0,
                'num_median_rel_error': float(statistics.median(rel_err_ov)) if len(rel_err_ov) > 0 else 0.0,
                'max_rel_error': float(max(rel_err_ov)) if len(rel_err_ov) > 0 else 0.0,
                'num_mean_abs_error_norm': float(statistics.mean(abs_err_ov)) / (float(max(map(abs, dg_ov_flat[idxs]))) + 1e-12) if len(abs_err_ov) > 0 else 0.0,
            },
            'G_hg': {
                'analytic_max_abs': float(max(map(abs, dg_hg_flat[idxs]))),
                'num_mean_abs_error': float(statistics.mean(abs_err_hg)) if len(abs_err_hg) > 0 else 0.0,
                'num_median_rel_error': float(statistics.median(rel_err_hg)) if len(rel_err_hg) > 0 else 0.0,
                'max_rel_error': float(max(rel_err_hg)) if len(rel_err_hg) > 0 else 0.0,
                'num_mean_abs_error_norm': float(statistics.mean(abs_err_hg)) / (float(max(map(abs, dg_hg_flat[idxs]))) + 1e-12) if len(abs_err_hg) > 0 else 0.0,
            }
        }

        print("有限差分校验结果:")
        print(f"  检查变量数: {result['n_checks']}")
        print(f"  G_ov 平均绝对误差: {result['G_ov']['num_mean_abs_error']:.6e}, 中位相对误差: {result['G_ov']['num_median_rel_error']:.6e}, 最大相对误差: {result['G_ov']['max_rel_error']:.6e}")
        print(f"  G_hg 平均绝对误差: {result['G_hg']['num_mean_abs_error']:.6e}, 中位相对误差: {result['G_hg']['num_median_rel_error']:.6e}, 最大相对误差: {result['G_hg']['max_rel_error']:.6e}")

        # 生成诊断信息：按绝对误差排序，输出 top-10（分别针对 G_ov 与 G_hg）
        import json as _json
        per_index_info_sorted_ov = sorted(per_index_info, key=lambda x: x['abs_err_ov'], reverse=True)
        per_index_info_sorted_hg = sorted(per_index_info, key=lambda x: x['abs_err_hg'], reverse=True)

        topk = 10
        diagnostics = {
            'summary': result,
            'top_errors_G_ov': per_index_info_sorted_ov[:topk],
            'top_errors_G_hg': per_index_info_sorted_hg[:topk]
        }

        # 保存诊断文件
        try:
            with open('fd_check_diagnostics.json', 'w') as _f:
                _json.dump(diagnostics, _f, indent=2)
            print(f"诊断信息已保存: fd_check_diagnostics.json (top {topk} 错误条目)")
        except Exception as _e:
            print(f"无法保存诊断信息: {_e}")

        return result
    
    def save_results(self, density, filename_base):
        """保存结果"""
        # 应用与优化过程中相同的滤波和投影
        # 注意：density 是设计变量 x，需要转换为物理密度 rho
        if isinstance(density, np.ndarray):
            x_tensor = torch.tensor(density, dtype=torch.float32, device=self.device)
        else:
            x_tensor = density
            
        with torch.no_grad():
            # 确保滤波器已设置（针对当前 density 的形状）
            # 如果使用了 ROI，这里的 density 是裁剪后的，shape 应该匹配 self.shape
            # 如果没匹配上，重新设置滤波器
            current_shape = tuple(x_tensor.shape)
            if not hasattr(self, 'filter_kernel') or (hasattr(self, 'shape') and tuple(self.shape) != current_shape):
                 self.setup_filter_gpu(current_shape)
            
            rho_final = self.filter_density_gpu(x_tensor)
            density_filtered = rho_final.cpu().numpy()

        # 保存密度场
        # 若有 ROI，回填到原始尺寸保存
        if hasattr(self, 'roi_slices') and self.roi_slices is not None and hasattr(self, 'original_shape') and self.original_shape is not None:
            full = np.zeros(self.original_shape, dtype=np.float32)
            sx, sy, sz = self.roi_slices
            full[sx, sy, sz] = density_filtered
            np.save(f"{filename_base}.npy", full)
            density_for_stl = full
        else:
            np.save(f"{filename_base}.npy", density_filtered)
            density_for_stl = density_filtered
        print(f"密度场已保存: {filename_base}.npy")
        
        # 使用连续密度生成 STL（更稳定、与优化一致）
        try:
            # STL 导出：若有 ROI 则用回填后的全局密度、无 origin 偏移
            self._export_stl_from_density(density_for_stl, f"{filename_base}.stl", level=0.5, step_size=1)
            print(f"优化结果保存完成")
        except Exception as e:
            print(f"❌ STL生成失败: {e}")

def main():
    """主函数"""
    #读参数
    parser = argparse.ArgumentParser(description='GPU拓扑优化器')
    
    # 接收输入文件路径 (绝对路径)
    parser.add_argument('--input', type=str, required=True, 
                        help='Input voxel(.npy) path')
    
    # 接收输出文件前缀 (绝对路径，例如 D:/res/result)
    parser.add_argument('--output', type=str, required=True, 
                        help='Output voxel(.npy) path')
    
    args = parser.parse_args()

    print("GPU加速拓扑优化器")
    print("="*60)
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")

    if not os.path.exists(args.input):
        print(f"错误: 找不到输入文件 -> {args.input}")
        sys.exit(1)
           
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        print(f"提示: 输出目录不存在，正在创建 -> {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

    print(f"output_dir: {output_dir}")

    if not torch.cuda.is_available():
        print("警告: CUDA不可用，将使用CPU")
        device = 'cpu'
    else:
        print(f"GPU加速可用: {torch.cuda.get_device_name()}")
        print(f"GPU内存: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
        device = 'cuda'
    
    # 加载输入
    #input_file = "voxelized_model.npy"
    #voxels = np.load(input_file)
    voxels = np.load(args.input)
    print(f"加载模型: {voxels.shape}")
    
    # 不降采样，使用原始精度
    original_shape = voxels.shape
    total_voxels = np.prod(original_shape)
    gpu_memory = torch.cuda.get_device_properties(0).total_memory / 1e9 if device == 'cuda' else 4.0
    
    print(f"原始模型: {original_shape}")
    print(f"总体素数: {total_voxels:,}")
    print(f"估算内存需求: {total_voxels * 4 * 8 / 1e9:.2f} GB")  # float32 * 8个张量
    
    if device == 'cuda' and total_voxels * 4 * 8 / 1e9 > gpu_memory * 0.8:
        print(f"警告: 内存需求可能超过GPU容量，建议使用分块处理")
        # 可以选择分块处理或者提示用户
        use_chunking = True
    else:
        use_chunking = False
    
    # 创建目标
    # 保持原始连续值，避免二值化造成的信息丢失
    target = voxels.astype(np.float32)
    
    # 初始化
    vol_frac = np.mean(target)
    initial_x = np.full(voxels.shape, vol_frac, dtype=np.float32)
    
    print(f"目标体积分数: {vol_frac:.4f}")
    
    # GPU优化
    start_time = time.time()
    optimizer = GPUTopologyOptimizer(
        print_direction=[1, 0, 0], 
        device=device, 
        use_chunking=use_chunking
    )
    

    #加默认迭代次数，确保硬约束下有足够收敛时间
    max_iter = 500

    # 默认不使用ROI裁剪，保持与baseline坐标对齐
    optimized_x, history = optimizer.optimize(initial_x, target, output_dir, max_iter=max_iter, use_roi=True)
    total_time = time.time() - start_time
    
    # 保存结果
    optimizer.save_results(optimized_x, output_dir + "/gpu_topology_optimized")
    
    print(f"\nGPU拓扑优化完成!")
    print(f"总时间: {total_time:.2f}s")
    print(f"平均每迭代: {total_time/len(history):.2f}s")
    print(f"最终目标函数: {history[-1]['objective']:.6f}")
    print(f"最终体积分数: {history[-1]['volume_fraction']:.4f}")
    # ... 在 optimizer.optimize(...) 调用之后 ...

    print("\n--- 开始有限差分校验 ---")
# 使用优化开始时的状态或优化后的状态进行检查
# 这里使用初始状态
    check_x = initial_x
    check_target = target

# 如果使用了ROI，需要用裁剪后的数据进行校验
    if optimizer.roi_slices:
        sx, sy, sz = optimizer.roi_slices
        check_x = initial_x[sx, sy, sz]
        check_target = target[sx, sy, sz]

    optimizer.finite_difference_check(check_x, check_target, n_checks=100)
    print("--- 有限差分校验结束 ---\n")

if __name__ == "__main__":
    main()
