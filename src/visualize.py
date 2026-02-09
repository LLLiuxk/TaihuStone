import numpy as np
import pyvista as pv
from scipy.ndimage import binary_erosion
from scipy.ndimage import label as cc_label

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
    voxel = np.load("D:/VSprojects/TaihuStone/origin_model.npy")
    voxel2 = np.load("D:/VSprojects/TaihuStone/opt_model.npy")
    voxel  = voxel  >= 0.5
    voxel2  = voxel2  >= 0.5
    visualize_A_reference(voxel, voxel2)
    visualize_B_reference(voxel, voxel2)
    #visualize_single_voxel(voxel)
    #visualize_voxel_diff(voxel, voxel2)
