import numpy as np
import pyvista as pv

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

def visualize_A_reference(A, B):
    A = A.astype(bool)
    B = B.astype(bool)
    assert A.shape == B.shape

    label = classify_voxels(A, B)

    grid = pv.ImageData()
    grid.dimensions = np.array(A.shape) + 1
    grid.spacing = (1, 1, 1)
    grid.origin = (0, 0, 0)

    grid.cell_data["label"] = label.flatten(order="F")

    plotter = pv.Plotter(title="A reference")

    both   = grid.extract_cells(grid.cell_data["label"] == 1)
    a_only = grid.extract_cells(grid.cell_data["label"] == 2)

    if both.n_cells > 0:
        plotter.add_mesh(both, color="lightgray", opacity=0.3)
    if a_only.n_cells > 0:
        plotter.add_mesh(a_only, color="red", opacity=1.0)

    plotter.show()
  
def visualize_B_reference(A, B):
    A = A.astype(bool)
    B = B.astype(bool)
    assert A.shape == B.shape

    label = classify_voxels(A, B)

    grid = pv.ImageData()
    grid.dimensions = np.array(A.shape) + 1
    grid.spacing = (1, 1, 1)
    grid.origin = (0, 0, 0)

    grid.cell_data["label"] = label.flatten(order="F")

    plotter = pv.Plotter(title="B reference")

    both   = grid.extract_cells(grid.cell_data["label"] == 1)
    b_only = grid.extract_cells(grid.cell_data["label"] == 3)

    if both.n_cells > 0:
        plotter.add_mesh(both, color="lightgray", opacity=0.3)
    if b_only.n_cells > 0:
        plotter.add_mesh(b_only, color="green", opacity=1.0)

    plotter.show()


if __name__ == "__main__":
    voxel = np.load("D:/VSprojects/TaihuStone/src/opt_model.npy")
    voxel2 = np.load("D:/VSprojects/TaihuStone/src/origin_model3.npy")
    voxel  = voxel  >= 0.5
    #visualize_A_reference(voxel, voxel2)
    #visualize_B_reference(voxel, voxel2)
    visualize_single_voxel(voxel)
    #visualize_voxel_diff(voxel, voxel2)
