import math
from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.geometry import Geometry
from adze_modeler.material import Material
from adze_modeler.metadata import Agros2DMetadata
from adze_modeler.metadata import FemmMetadata
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from copy import copy
from pathlib import Path
from random import uniform
from shutil import rmtree
from time import ctime
from uuid import uuid4

from numpy import linspace
from numpy import meshgrid

export_location = Path(__file__).parent / "snapshots"


def build(platform, X, msize, cleanup=True, customid=None):

    p = copy(platform)
    model_id = customid or str(uuid4())

    print(ctime(), model_id, end="....BUILDING...")
    modelpath = export_location / model_id
    if not export_location.exists():
        export_location.mkdir(parents=True, exist_ok=True)

    modelpath.mkdir(parents=False, exist_ok=True)
    with open(modelpath / "X.csv", "w") as f:
        f.write(",".join([str(i) for i in X]))

    p.metadata.file_script_name = str(modelpath / p.metadata.file_script_name)
    p.metadata.file_metrics_name = str(modelpath / p.metadata.file_metrics_name)

    mp_bound = ModelPiece("bound")
    mp_bound.load_piece_from_svg("problem_boundary.svg")

    mp_control = ModelPiece("core")
    mp_control.load_piece_from_svg("core.svg")
    mp_control.put(0, -5)

    mp_coil = ModelPiece("coil")
    mp_coil.load_piece_from_svg("coil.svg")

    snapshot = Snapshot(p)
    exctitation = Material("J+")
    exctitation.Je = 2e6
    air = Material("air")
    control = Material("control")

    exctitation.meshsize = msize
    air.meshsize = msize
    control.meshsize = msize

    snapshot.add_material(exctitation)
    snapshot.add_material(air)
    snapshot.add_material(control)

    snapshot.assign_material(3, 0, name="control")
    snapshot.assign_material(30, 30, name="air")

    b1 = DirichletBoundaryCondition(
        name="a0", field_type=snapshot.platform.metadata.problem_type, magnetic_potential=0.0
    )

    N = len(X)
    h = 1.5
    geom = Geometry()
    geom.merge_geometry(mp_bound.geom)
    geom.merge_geometry(mp_control.geom)

    for i, ri in enumerate(X):
        coil = mp_coil.spawn()
        offsetz = i * h - N * h / 2
        coil.put(ri, offsetz)
        geom.merge_geometry(coil.geom)
        snapshot.assign_material(ri + h / 2, offsetz + h / 2, name="J+")

    geom.generate_intersections()
    snapshot.add_geometry(geom)
    snapshot.add_boundary_condition(b1)

    snapshot.assign_boundary_condition(0, 0, "a0")
    snapshot.assign_boundary_condition(0, -20, "a0")
    snapshot.assign_boundary_condition(0, 20, "a0")
    snapshot.assign_boundary_condition(35, 35, "a0")
    snapshot.assign_boundary_condition(35, -35, "a0")
    snapshot.assign_boundary_condition(70, 0, "a0")

    geom.export_geom(modelpath / f"{model_id}.svg")

    Nx = 10
    Ny = 10
    px = linspace(0, 5, Nx)
    py = linspace(-5, 5, Ny)
    xv, yv = meshgrid(px, py, sparse=False, indexing="xy")

    for i in range(Nx):
        for j in range(Ny):
            eval_point = (xv[j, i], yv[j, i])
            snapshot.add_postprocessing("point_value", eval_point, "Bz")
            snapshot.add_postprocessing("point_value", eval_point, "Br")

    snapshot.add_postprocessing("mesh_info", None, None)

    print("EXPORTING", end="...")
    snapshot.export()
    print(f"EXECUTING({snapshot.platform.metadata.compatible_platform})", end="...")
    try:
        snapshot.execute()
    except Exception as e:
        print("FAILED TO SOLVE")
        return None
    else:
        print("SOLVED", end="...EXTRACTING SOL.")
    try:
        res = snapshot.retrive_results()
        res["platform"] = snapshot.platform.metadata.compatible_platform
        print("DONE")
    except FileNotFoundError:
        print("DISCARDED (Invalid config)")
        return None

    if cleanup:
        rmtree(modelpath)

    return res
