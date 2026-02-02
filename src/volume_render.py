import numpy as np
import pyvista as pv


def render_sdf_volume(
    npy_path,
    spacing=(1.0, 1.0, 1.0),
    clim=None,
):
    """
    Volume rendering for SDF stored as .npy

    Parameters
    ----------
    npy_path : str
        Path to sdf.npy, shape = (nx, ny, nz)
    spacing : tuple(float, float, float)
        Voxel size (dx, dy, dz)
    clim : tuple(float, float) or None
        Color limit for SDF visualization
    """

    # 1. load sdf
    sdf = np.load(npy_path).astype(np.float32)
    nx, ny, nz = sdf.shape

    print(f"[INFO] SDF shape: {sdf.shape}")
    print(f"[INFO] SDF min/max: {sdf.min():.4f}, {sdf.max():.4f}")

    # 2. create VTK ImageData
    grid = pv.ImageData(
        dimensions=(nx, ny, nz),
        spacing=spacing,
        origin=(-0.5, -0.5, -0.5),
    )

    # VTK expects Fortran order
    grid["sdf"] = sdf.flatten(order="F")

    # 3. plotter
    pl = pv.Plotter()
    pl.set_background("white")

    opacity = [0.0, 0.1, 0.3, 0.6, 1.0]
    # [
    #     (-0.2, 0.0),
    #     (-0.05, 0.1),
    #     (0.0, 0.9),
    #     (0.05, 0.1),
    #     (0.7, 0.0)
    # ]

    pl.add_volume(
        grid,
        scalars="sdf",
        cmap="coolwarm",
        opacity="sigmoid",
        clim= clim,
        shade=True,

    )
    #
    # # 4. volume rendering
    # pl.add_volume(
    #     grid,
    #     scalars="sdf",
    #     cmap="coolwarm",
    #     clim=clim,
    #     opacity="sigmoid",   # very suitable for SDF
    #     shade=True,
    # )

    # optional: add axes
    pl.add_axes()
    pl.show()


if __name__ == "__main__":
    render_sdf_volume(
        "origin_model.npy",
        spacing=(1.0, 1.0, 1.0),
        clim=(-0.5, 0.5),
    )
