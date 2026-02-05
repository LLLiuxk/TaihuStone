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
        color="lightgray",
        show_edges=False,
        opacity=0.1
    )
    plotter.show()


if __name__ == "__main__":
    voxel = np.load("D:/VSprojects/TaihuStone/src/origin_model.npy")
    visualize_single_voxel(voxel)
