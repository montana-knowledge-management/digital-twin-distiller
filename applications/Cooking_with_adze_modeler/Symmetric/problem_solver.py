import math
import operator
from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.boundaries import NeumannBoundaryCondition
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

current_dir = Path(__file__).parent
export_location = Path(__file__).parent / "snapshots"


def build(platform, X, cleanup=True, customid=None):
    p = copy(platform)
    model_id = customid or str(uuid4())

    print(ctime(), model_id, end="....BUILDING...")
    modelpath = export_location / model_id
    if not export_location.exists():
        export_location.mkdir(parents=True, exist_ok=True)

    modelpath.mkdir(parents=True, exist_ok=True)
    with open(modelpath / "X.csv", "w") as f:
        f.write(",".join([str(i) for i in X]))

    p.metadata.file_script_name = str(modelpath / p.metadata.file_script_name)
    p.metadata.file_metrics_name = str(modelpath / p.metadata.file_metrics_name)

    mp_bound = ModelPiece("bound")
    mp_bound.load_piece_from_svg(current_dir / "problem_boundary.svg")
    mp_bound.put(0, 0)

    mp_control = ModelPiece("core")
    mp_control.load_piece_from_svg(current_dir / "core.svg")
    mp_control.put(0, 0)

    mp_coil = ModelPiece("coil")
    mp_coil.load_piece_from_svg(current_dir / "coil.svg")

    snapshot = Snapshot(p)
    exctitation = Material("J+")
    exctitation.Je = 2e6
    air = Material("air")
    control = Material("control")

    # exctitation.meshsize = 1
    air.meshsize = 0.5
    control.meshsize = 0.1

    snapshot.add_material(exctitation)
    snapshot.add_material(air)
    snapshot.add_material(control)

    snapshot.assign_material(3, 1, name="control")
    snapshot.assign_material(30, 30, name="air")

    b1 = DirichletBoundaryCondition(
        name="a0", field_type=snapshot.platform.metadata.problem_type, magnetic_potential=0.0
    )
    n0 = NeumannBoundaryCondition(name="n0", field_type=snapshot.platform.metadata.problem_type, surface_current=0.0)

    N = len(X)
    h = 1.5
    geom = Geometry()
    geom.merge_geometry(mp_bound.geom)
    geom.merge_geometry(mp_control.geom)

    for i, ri in enumerate(X):
        coil = mp_coil.spawn()
        offsetz = i * h
        coil.put(ri, offsetz)
        geom.merge_geometry(coil.geom)
        snapshot.assign_material(ri + h / 2, offsetz + h / 2, name="J+")

    geom.generate_intersections()
    snapshot.add_geometry(geom)
    snapshot.add_boundary_condition(b1)
    snapshot.add_boundary_condition(n0)

    snapshot.assign_boundary_condition(0, 2.5, "a0")
    snapshot.assign_boundary_condition(0, 25, "a0")
    snapshot.assign_boundary_condition(40, 140, "a0")
    snapshot.assign_boundary_condition(140, 40, "a0")

    snapshot.assign_boundary_condition(2, 0, "n0")
    snapshot.assign_boundary_condition(5.001, 0, "n0")
    snapshot.assign_boundary_condition(130, 0, "n0")
    snapshot.assign_boundary_condition(X[0] + 0.5, 0, "n0")

    geom.export_geom(modelpath / f"{model_id}.svg")

    Nx = 11
    Ny = 6
    px = linspace(0.001, 4.9, Nx)
    py = linspace(0.001, 4.9, Ny)
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
        snapshot.execute(timeout=1000)
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


def calculate_objective_functions(X, resn, res, resp):
    Bz = [pointvalue[2] for pointvalue in res["Bz"]]  # [x, y, Bz(x, y)]
    Br = [pointvalue[2] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
    xi = [pointvalue[0] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
    yi = [pointvalue[1] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
    nb_nodes = res["nodes"]

    # Calculate F1
    B0 = 2e-3
    F1 = max(map(lambda Bz_i: abs(Bz_i - B0), Bz))

    # Calculate F2

    Bzp = [pointvalue[2] for pointvalue in resp["Bz"]]
    Brp = [pointvalue[2] for pointvalue in resp["Br"]]

    Bzn = [pointvalue[2] for pointvalue in resn["Bz"]]
    Brn = [pointvalue[2] for pointvalue in resn["Br"]]

    deltaBpz = map(operator.abs, map(operator.sub, Bzp, Bz))
    deltaBpr = map(operator.abs, map(operator.sub, Brp, Br))
    deltaBp = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBpz, deltaBpr))

    deltaBnz = map(operator.abs, map(operator.sub, Bzn, Bz))
    deltaBnr = map(operator.abs, map(operator.sub, Brn, Br))
    deltaBn = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBnz, deltaBnr))

    F2 = max(map(operator.add, deltaBp, deltaBn))

    # Calcukate F3
    F3 = sum(X)

    return [F1, F2, F3]


if __name__ == "__main__":
    from random import seed, uniform

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
    agros_metadata.nb_refinements = 0
    agros_metadata.polyorder = 2
    agros_metadata.adaptivity = "hp-adaptivity"
    agros_metadata.adaptivity_tol = 0.001

    platform_femm = Femm(femm_metadata)
    platform_agros = Agros2D(agros_metadata)

    # seed(42)
    # X = [uniform(5.5, 50) for _ in range(4)]
    # X.extend([uniform(1, 50) for _ in range(6)])
    #
    # X = [10] * 10
    # resf = build(platform_femm, X, cleanup=False)
    #
    # # resf = build(platform_agros, X, cleanup=False)
    # print(resf)

    # X = [10]*10
    # fd1, fd2, fd3 = 1.36e-4, 5.59e-5, 100
    # f1 = 9.40e-05, f2 = 6.77e-05, f3 = 1.00e+02
    # f1: -30.85189826773091
    # f2: 21.07567781921821
    # f3: 0.0

    X = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    fd1, fd2, fd3 = 8.18e-4, 3.01e-4, 105
    # f1 = 8.64e-04, f2 = 3.01e-04, f3 = 1.05e+02
    # f1: 5.570869222493864
    # f2: 0.1035562092624526
    # f3: 0.0

    # f1 = 8.36e-04, f2 = 2.93e-04, f3 = 1.05e+02
    # f1: 2.252706699266494
    # f2: -2.6387197953922583
    # f3: 0.0

    # X = [8.08, 14.9, 6.74, 16.7, 5.45, 10.6, 11.7, 11.1, 13.6, 7.99]
    # fd1, fd2, fd3 = 3.13e-5, 8.54e-5, 106.8
    # f1 = 6.80e-05, f2 = 8.89e-05, f3 = 1.07e+02
    # f1: 117.10884483035449
    # f2: 4.050306334556429
    # f3: 0.0561797752808877

    X = [7.71, 9.93, 16.9, 6.19, 18.3, 8.41, 7.19, 14.1, 8.91, 10.32]
    fd1, fd2, fd3 = 7.28e-5, 1.91e-4, 108
    # FEMM
    # f1 = 8.60e-05, f2 = 1.77e-04, f3 = 1.08e+02
    # f1: 18.085674698189678
    # f2: -7.571712770424163
    # f3: -0.037037037037055985

    # Agros2D
    # f1 = 8.36e-05, f2 = 1.78e-04, f3 = 1.08e+02
    # f1: 14.808593749999806
    # f2: -6.6372475683339305
    # f3: -0.037037037037055985

    # X = [11, 8, 14, 10, 6, 11, 8, 11, 13, 13]
    # fd1, fd2, fd3 = 3.16e-5, 1.21e-5, 105
    # f1 = 8.06e-05, f2 = 9.76e-05, f3 = 1.05e+02
    # f1: 155.06827567285455
    # f2: 706.7588240828788
    # f3: 0.0

    # X = [31.3, 29.7, 33.9, 33.0, 20.6, 16.3, 12.4, 5.53, 13.2, 5.54]
    # fd1, fd2, fd3 = 9.45e-4, 3.87e-6, 201.4
    # f1 = 9.04e-04, f2 = 6.06e-06, f3 = 2.01e+02
    # f1: -4.309629550657461
    # f2: 56.51081097301244
    # f3: 0.03475670307844746

    # X = [16.1, 16.2, 30.9, 14.3, 12.5, 11.5, 25.1, 5.5, 44.3, 6.96]
    # fd1, fd2, fd3 = 5.77e-4, 2.98e-5, 183.3
    # f1 = 5.43e-04, f2 = 2.73e-05, f3 = 1.83e+02
    # f1: -5.822371325167932
    # f2: -8.260566352371026
    # f3: 0.03273322422257166
    #
    # X = [29, 47, 45, 33, 36, 15, 29, 31, 6, 40]
    # fd1, fd2, fd3 = 1.08e-3, 1.52e-5, 311
    # f1 = 1.01e-03, f2 = 1.38e-05, f3 = 3.11e+02
    # f1: -6.107418505092587
    # f2: -9.117013622394852
    # f3: 0.0

    # X = []
    # fd1, fd2, fd3 =

    # X = list(reversed(X))
    Xp = [xi + 0.5 for xi in X]
    Xn = [xi - 0.5 for xi in X]
    res_m = build(platform_agros, Xn, cleanup=False, customid="alma_minus")
    resa_0 = build(platform_agros, X, cleanup=False, customid="alma_0")
    resa_p = build(platform_agros, Xp, cleanup=False, customid="alma_plus")

    f1, f2, f3 = calculate_objective_functions(X, res_m, resa_0, resa_p)

    print(f"{f1=:.2e}, {f2=:.2e}, {f3=:.2e}")
    print("f1:", (f1 - fd1) / fd1 * 100)
    print("f2:", (f2 - fd2) / fd2 * 100)
    print("f3:", (f3 - fd3) / fd3 * 100)
    """
    B0 = 2e-3
    print((f1-B0)/B0*100)
    """
