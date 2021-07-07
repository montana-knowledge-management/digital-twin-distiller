import math
from uuid import uuid4
from adze_modeler.metadata import Agros2DMetadata, FemmMetadata
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.material import Material
from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.geometry import Geometry
from pathlib import Path
from numpy import meshgrid, linspace
from shutil import rmtree
from copy import copy
from time import ctime
from random import uniform

current_dir = Path(__file__).parent
export_location = Path(__file__).parent / "snapshots"

def build(platform, X, cleanup=True, customid=None):

    p = copy(platform)
    model_id = customid or str(uuid4())

    print( ctime(), model_id, end='....BUILDING...')
    modelpath = export_location / model_id
    if not export_location.exists():
        export_location.mkdir(parents=True, exist_ok=True)

    modelpath.mkdir(parents=True, exist_ok=True)
    with open(modelpath / "X.csv", "w") as f:
        f.write(','.join([str(i) for i in X]))

    p.metadata.file_script_name = str(modelpath / p.metadata.file_script_name)
    p.metadata.file_metrics_name = str(modelpath / p.metadata.file_metrics_name)

    mp_bound = ModelPiece("bound")
    mp_bound.load_piece_from_svg(current_dir / 'problem_boundary.svg')

    mp_control = ModelPiece("core")
    mp_control.load_piece_from_svg(current_dir / 'core.svg')
    mp_control.put(0, -5)

    mp_coil = ModelPiece("coil")
    mp_coil.load_piece_from_svg(current_dir / 'coil.svg')


    snapshot = Snapshot(p)
    exctitation = Material("J+")
    exctitation.Je = 2e6
    air = Material("air")
    control = Material("control")

    exctitation.meshsize = 1
    air.meshsize = 0.7
    control.meshsize = 0.1

    snapshot.add_material(exctitation)
    snapshot.add_material(air)
    snapshot.add_material(control)

    snapshot.assign_material(3, 0, name="control")
    snapshot.assign_material(30, 30, name="air")

    b1 = DirichletBoundaryCondition(name="a0", field_type=snapshot.platform.metadata.problem_type, magnetic_potential=0.0)

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
        snapshot.assign_material(ri + h/2, offsetz + h/2, name="J+")

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
    px = linspace(0.001, 5, Nx)
    py = linspace(-5, 5, Ny)
    xv, yv = meshgrid(px, py, sparse=False, indexing='xy')

    for i in range(Nx):
        for j in range(Ny):
            eval_point = (xv[j, i], yv[j, i])
            snapshot.add_postprocessing('point_value', eval_point, 'Bz')
            snapshot.add_postprocessing('point_value', eval_point, 'Br')

    snapshot.add_postprocessing("mesh_info", None, None)

    print("EXPORTING", end='...')
    snapshot.export()
    print(f"EXECUTING({snapshot.platform.metadata.compatible_platform})", end='...')
    try:
        snapshot.execute()
    except Exception as e:
        print("FAILED TO SOLVE")
        return None
    else:
        print("SOLVED", end='...EXTRACTING SOL.')
    try:
        res = snapshot.retrive_results()
        res['platform'] = snapshot.platform.metadata.compatible_platform
        print("DONE")
    except FileNotFoundError:
        print("DISCARDED (Invalid config)")
        return None

    if cleanup:
        rmtree(modelpath)

    return res



if __name__=='__main__':
    femm_metadata = FemmMetadata()
    femm_metadata.problem_type = "magnetic"
    femm_metadata.coordinate_type = "axisymmetric"
    femm_metadata.file_script_name = "femm_solver_script"
    femm_metadata.file_metrics_name = "femm_solution.csv"
    femm_metadata.unit = "millimeters"
    femm_metadata.smartmesh = False

    agros_metadata = Agros2DMetadata()
    agros_metadata.file_script_name = "agros_solver_script"
    agros_metadata.file_metrics_name = "agros_solution.csv"
    agros_metadata.problem_type = "magnetic"
    agros_metadata.coordinate_type = "axisymmetric"
    agros_metadata.analysis_type = "steadystate"
    agros_metadata.unit = 1e-3
    agros_metadata.nb_refinements = 2
    agros_metadata.adaptivity = "hp-adaptivity"
    agros_metadata.adaptivity_tol = 0.2

    platform_femm = Femm(femm_metadata)
    platform_agros = Agros2D(agros_metadata)
    from random import seed, uniform
    seed(42)
    X = [uniform(1, 50) for _ in range(6)]
    X.extend([uniform(5, 50) for _ in range(12)])
    X.extend([uniform(1, 50) for _ in range(6)])

    for i in range(10):
        resf = build(platform_femm, X, cleanup=False, customid="nodes-F1")
        resa = build(platform_agros, X, cleanup=False, customid="nodes-F1")