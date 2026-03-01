import numpy as np
import pyvista as pv
from scipy.ndimage import binary_erosion
from scipy.ndimage import label as cc_label
from scipy.ndimage import maximum_filter

def compute_support_voxels(voxel):
    """
    voxel: shape (nx, ny, nz)
    打印方向 = +z = axis 2
    """
    voxel = voxel.astype(bool)
    nx, ny, nz = voxel.shape

    support_mask = np.zeros_like(voxel, dtype=bool)

    for x in range(1, nx):

        current_layer = voxel[x, :, :]      # y-z 平面
        below_layer   = voxel[x-1, :, :]

        # 在 y-z 平面做 3×3 邻域
        support_region = maximum_filter(below_layer, size=3)

        unsupported = current_layer & (~support_region)

        support_mask[x, :, :] = unsupported

    return support_mask


def compute_support_voxels2(voxel):
    """
    voxel: 3D bool array
    return:
        support_mask: 需要支撑的体素 (bool)
    """

    voxel = voxel.astype(bool)
    nx, ny, nz = voxel.shape

    support_mask = np.zeros_like(voxel, dtype=bool)

    # 从第 1 层开始（第 0 层默认认为有底板支撑）
    for k in range(1, nz):
        current_layer = voxel[:, :, k]
        below_layer   = voxel[:, :, k-1]

        # 对 below_layer 做 3x3 邻域扩张（卷积思想）
        padded = np.pad(below_layer, 1, mode='constant', constant_values=False)

        support_region = np.zeros_like(below_layer, dtype=bool)

        for dx in range(3):
            for dy in range(3):
                support_region |= padded[dx:dx+nx, dy:dy+ny]

        # 当前层存在实体 且 下层3x3无支撑
        unsupported = current_layer & (~support_region)

        support_mask[:, :, k] = unsupported

    return support_mask


def visualize_support(voxel):
    voxel = voxel.astype(bool)

    support_mask = compute_support_voxels(voxel)

    total_voxels = np.count_nonzero(voxel)
    support_voxels = np.count_nonzero(support_mask)

    print("总实体体素数:", total_voxels)
    print("需要支撑体素数:", support_voxels)
    print("支撑比例:", support_voxels / total_voxels)

    # 构建 grid
    grid = pv.ImageData()
    grid.dimensions = np.array(voxel.shape) + 1
    grid.spacing = (1, 1, 1)

    label = np.zeros_like(voxel, dtype=np.uint8)
    label[voxel] = 1
    label[support_mask] = 2

    grid.cell_data["label"] = label.flatten(order="F")

    plotter = pv.Plotter()

    normal_voxels  = grid.extract_cells(grid.cell_data["label"] == 1)
    support_voxels = grid.extract_cells(grid.cell_data["label"] == 2)

    if normal_voxels.n_cells > 0:
        plotter.add_mesh(normal_voxels, color="gray", opacity=0.3)

    if support_voxels.n_cells > 0:
        plotter.add_mesh(support_voxels, color="red", opacity=1.0)

    plotter.show()

def visualize_single_voxel(voxel, threshold=0):
    """
    voxel: (nx, ny, nz) numpy array
    """
    occ = voxel > threshold

    nx, ny, nz = occ.shape

    # 生成 StructuredGrid
    grid = pv.ImageData()
    grid.dimensions = np.array(occ.shape) + 1
    grid.spacing = (1.0, 1.0, 1.0)
    grid.origin = (0, 0, 0)

    # Cell data
    grid.cell_data["occ"] = occ.flatten(order="F")

    # 提取 occupied cells
    solid = grid.extract_cells(grid.cell_data["occ"])

    plotter = pv.Plotter()
    plotter.add_mesh(
        solid,
        color="gray",
        show_edges=False,
        opacity=0.8
    )
    plotter.show()

def classify_voxels(A, B):
    """
    A, B: boolean occupancy grids
    return: label grid
    0 = empty
    1 = A ∩ B
    2 = A - B
    3 = B - A
    """
    label = np.zeros_like(A, dtype=np.uint8)
    label[(A == 1) & (B == 1)] = 1
    label[(A == 1) & (B == 0)] = 2
    label[(A == 0) & (B == 1)] = 3
    return label
    

def visualize_voxel_diff(A, B):
    A = A.astype(bool)
    B = B.astype(bool)
    assert A.shape == B.shape

    label = classify_voxels(A, B)

    grid = pv.ImageData()
    grid.dimensions = np.array(A.shape) + 1
    grid.spacing = (1, 1, 1)
    grid.origin = (0, 0, 0)

    grid.cell_data["label"] = label.flatten(order="F")

    plotter = pv.Plotter()

    both   = grid.extract_cells(grid.cell_data["label"] == 1)
    a_only = grid.extract_cells(grid.cell_data["label"] == 2)
    b_only = grid.extract_cells(grid.cell_data["label"] == 3)

    if both.n_cells > 0:
        plotter.add_mesh(both, color="lightgray", opacity=0.3)
    if a_only.n_cells > 0:
        plotter.add_mesh(a_only, color="red", opacity=0.50)
    if b_only.n_cells > 0:
        plotter.add_mesh(b_only, color="blue", opacity=0.50)

    plotter.show()

def filter_small_components(mask, min_size):
    labeled, num = cc_label(mask)
    sizes = np.bincount(labeled.ravel())

    keep = sizes >= min_size
    keep[0] = False  # background

    return keep[labeled]

def classify_voxels_connected(A, B, min_size=50):
    label = np.zeros_like(A, dtype=np.uint8)

    both   = (A == 1) & (B == 1)
    a_only = (A == 1) & (B == 0)
    b_only = (A == 0) & (B == 1)

    a_only_big = filter_small_components(a_only, min_size)
    b_only_big = filter_small_components(b_only, min_size)

    label[both] = 1
    label[a_only_big] = 2
    label[b_only_big] = 3

    return label

def surface_voxels(mask):
    """返回 mask 的表面体素"""
    eroded = binary_erosion(mask)
    return mask & (~eroded)

def classify_voxels_surface(A, B):
    label = np.zeros_like(A, dtype=np.uint8)

    # 原始三类
    both   = (A == 1) & (B == 1)
    a_only = (A == 1) & (B == 0)
    b_only = (A == 0) & (B == 1)

    # 只保留“表面相关”的变化
    A_surface = surface_voxels(A)
    B_surface = surface_voxels(B)

    a_only_surface = a_only & A_surface
    b_only_surface = b_only & B_surface

    label[both] = 1
    label[a_only_surface] = 2
    label[b_only_surface] = 3

    return label



def visualize_A_reference(A, B):
    A = A.astype(bool)
    B = B.astype(bool)
    assert A.shape == B.shape

    #label = classify_voxels(A, B)
    #label = classify_voxels_surface(A, B)
    label = classify_voxels_connected(A, B, min_size=3)

    grid = pv.ImageData()
    grid.dimensions = np.array(A.shape) + 1
    grid.spacing = (1, 1, 1)
    grid.origin = (0, 0, 0)

    grid.cell_data["label"] = label.flatten(order="F")

    plotter = pv.Plotter(title="A reference")

    both   = grid.extract_cells(grid.cell_data["label"] == 1)
    a_only = grid.extract_cells(grid.cell_data["label"] == 2)

    if both.n_cells > 0:
        plotter.add_mesh(both, color="gray", opacity=0.4)
    if a_only.n_cells > 0:
        plotter.add_mesh(a_only, color="red", opacity=0.80)

    plotter.show()
  
def visualize_B_reference(A, B):
    A = A.astype(bool)
    B = B.astype(bool)
    assert A.shape == B.shape

    #label = classify_voxels(A, B)
    #label = classify_voxels_surface(A, B)
    label = classify_voxels_connected(A, B, min_size=15)

    grid = pv.ImageData()
    grid.dimensions = np.array(A.shape) + 1
    grid.spacing = (1, 1, 1)
    grid.origin = (0, 0, 0)

    grid.cell_data["label"] = label.flatten(order="F")

    plotter = pv.Plotter(title="B reference")

    both   = grid.extract_cells(grid.cell_data["label"] == 1)
    b_only = grid.extract_cells(grid.cell_data["label"] == 3)

    if both.n_cells > 0:
        plotter.add_mesh(both, color="gray", opacity=0.4)
    if b_only.n_cells > 0:
        plotter.add_mesh(b_only, color="green", opacity=0.80)

    plotter.show()


if __name__ == "__main__":
    # voxel = np.load("D:/VSprojects/TaihuStone/origin_model.npy")
    # voxel2 = np.load("D:/VSprojects/TaihuStone/opt_model.npy")
    # voxel  = voxel  >= 0.5
    # voxel2  = voxel2  >= 0.5
    # visualize_A_reference(voxel, voxel2)
    # visualize_B_reference(voxel, voxel2)

    rho = np.load("D:/VSprojects/TaihuStone/test.npy")
    voxel = rho >= 0.5
    visualize_support(voxel)
    #visualize_single_voxel(voxel)
    #visualize_voxel_diff(voxel, voxel2)
