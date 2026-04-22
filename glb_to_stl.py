from __future__ import annotations

import sys
from pathlib import Path


def require_trimesh():
    try:
        import trimesh  # type: ignore
    except ImportError:
        print("Error: Python package 'trimesh' is required.")
        print("Install it with:")
        print(f'  "{sys.executable}" -m pip install trimesh numpy')
        raise SystemExit(1)
    return trimesh


def collect_inputs(arguments: list[str]) -> list[Path]:
    inputs: list[Path] = []

    if not arguments:
        arguments = [str(Path.cwd())]

    for item in arguments:
        path = Path(item).expanduser()
        if path.is_dir():
            inputs.extend(sorted(path.glob("*.glb")))
        elif path.is_file() and path.suffix.lower() == ".glb":
            inputs.append(path)
        else:
            print(f"Skip: {path} is not a .glb file or directory.")

    unique: list[Path] = []
    seen: set[Path] = set()
    for path in inputs:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(resolved)
    return unique


def scene_to_mesh(trimesh, loaded):
    if isinstance(loaded, trimesh.Trimesh):
        return loaded

    if not isinstance(loaded, trimesh.Scene):
        raise ValueError(f"Unsupported GLB content type: {type(loaded).__name__}")

    meshes = []
    for geometry in loaded.dump():
        if isinstance(geometry, trimesh.Trimesh) and len(geometry.faces) > 0:
            meshes.append(geometry)

    if not meshes:
        raise ValueError("No triangle mesh geometry was found.")

    if len(meshes) == 1:
        return meshes[0]

    return trimesh.util.concatenate(meshes)


def convert_file(trimesh, glb_path: Path) -> Path:
    stl_path = glb_path.with_suffix(".stl")
    loaded = trimesh.load(glb_path, force="scene", process=False)
    mesh = scene_to_mesh(trimesh, loaded)

    if len(mesh.vertices) == 0 or len(mesh.faces) == 0:
        raise ValueError("Mesh has no vertices or faces.")

    mesh.export(stl_path, file_type="stl")
    return stl_path


def main() -> int:
    trimesh = require_trimesh()
    inputs = collect_inputs(sys.argv[1:])

    if not inputs:
        print("No .glb files found.")
        return 1

    failed = 0
    for glb_path in inputs:
        try:
            output = convert_file(trimesh, glb_path)
            print(f"OK: {glb_path} -> {output}")
        except Exception as exc:
            failed += 1
            print(f"FAIL: {glb_path}")
            print(f"      {exc}")

    print()
    print(f"Done. Success: {len(inputs) - failed}, Failed: {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
