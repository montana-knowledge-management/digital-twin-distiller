import functools
import multiprocessing
from adze_modeler.boundaries import AntiPeriodicBoundaryCondition
from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.geometry import Geometry
from adze_modeler.material import Material
from adze_modeler.metadata import FemmMetadata
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.objects import CubicBezier
from adze_modeler.objects import Line
from adze_modeler.objects import Node
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from copy import copy
from pathlib import Path
from time import perf_counter

from numpy import linspace

basepath = Path(__file__).parent
exportpath = basepath / "snapshots"
exportpath.mkdir(exist_ok=True)


def generate_snapshot(stator, rotor, X, di, export_loc=None):
    export_loc = export_loc or exportpath
    femm_metadata = FemmMetadata()
    femm_metadata.problem_type = "magnetic"
    femm_metadata.coordinate_type = "planar"
    femm_metadata.file_script_name = export_loc / f"femm_solver_script_{int(di * 1000)}"
    femm_metadata.file_metrics_name = export_loc / f"femm_solution_{int(di * 1000)}.csv"
    femm_metadata.unit = "millimeters"
    femm_metadata.smartmesh = True
    femm_metadata.depth = 200
    platform = Femm(femm_metadata)

    snapshot = Snapshot(platform)

    a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)
    pb1 = AntiPeriodicBoundaryCondition("PB1", field_type="magnetic")
    pb2 = AntiPeriodicBoundaryCondition("PB2", field_type="magnetic")
    pb3 = AntiPeriodicBoundaryCondition("PB3", field_type="magnetic")
    pb4 = AntiPeriodicBoundaryCondition("PB4", field_type="magnetic")
    pb5 = AntiPeriodicBoundaryCondition("PB5", field_type="magnetic")
    pb6 = AntiPeriodicBoundaryCondition("PB6", field_type="magnetic")
    pb7 = AntiPeriodicBoundaryCondition("PB7", field_type="magnetic")

    snapshot.add_boundary_condition(a0)
    snapshot.add_boundary_condition(pb1)
    snapshot.add_boundary_condition(pb2)
    snapshot.add_boundary_condition(pb3)
    snapshot.add_boundary_condition(pb4)
    snapshot.add_boundary_condition(pb5)
    snapshot.add_boundary_condition(pb6)
    snapshot.add_boundary_condition(pb7)

    magnet = Material("N38")
    # magnet.mu_r = 1.05
    magnet.coercivity = 944771
    magnet.conductivity = 0.667 * 1e6
    magnet.remanence_angle = 90

    air = Material("air")

    steel = Material("M-27")
    steel.b = [
        0,
        0.05,
        0.1,
        0.15,
        0.2,
        0.25,
        0.3,
        0.35,
        0.4,
        0.45,
        0.5,
        0.55,
        0.6,
        0.65,
        0.7,
        0.75,
        0.8,
        0.85,
        0.9,
        0.95,
        1,
        1.05,
        1.1,
        1.15,
        1.2,
        1.25,
        1.3,
        1.35,
        1.4,
        1.45,
        1.5,
        1.55,
        1.6,
        1.65,
        1.7,
        1.75,
        1.8,
        1.85,
        1.9,
        1.95,
        2,
        2.05,
        2.1,
        2.15,
        2.2,
        2.25,
        2.3,
    ]
    steel.h = [
        0,
        17.059828,
        25.634971,
        31.338354,
        35.778997,
        39.602124,
        43.123143,
        46.520439,
        49.908177,
        53.368284,
        56.966318,
        60.760271,
        64.806105,
        69.161803,
        73.890922,
        79.066315,
        84.774676,
        91.122722,
        98.246299,
        106.324591,
        115.603417,
        126.435124,
        139.349759,
        155.187082,
        175.350538,
        202.312017,
        240.640455,
        299.118027,
        394.993386,
        561.726177,
        859.328763,
        1375.466888,
        2191.246914,
        3328.145908,
        4760.506172,
        6535.339449,
        8788.970657,
        11670.804347,
        15385.186211,
        20246.553031,
        26995.131141,
        38724.496369,
        64917.284463,
        101489.309338,
        137202.828961,
        176835.706764,
        216374.283609,
    ]
    steel.conductivity = 2 * 1e6
    steel.thickness = 0.635
    steel.lamination_type = "inplane"
    steel.fill_factor = 0.98

    snapshot.add_material(magnet)
    snapshot.add_material(air)
    snapshot.add_material(steel)

    snapshot.assign_material(32.8, 82.85, name="M-27")
    snapshot.assign_material(di + 20, -30, name="M-27")

    snapshot.assign_material(3, 37, name="air")
    snapshot.assign_material(23.6, 36.3, name="air")
    snapshot.assign_material(44.7, 36.8, name="air")
    snapshot.assign_material(65.1, 36.9, name="air")
    snapshot.assign_material(0.08, 0.69, name="air")
    snapshot.assign_material(di + 3, -6, name="air")
    snapshot.assign_material(di + 65, -6, name="air")
    snapshot.assign_material(di + 32, -6, name="N38")

    rotor.put(di, -50.75)
    geom = Geometry()
    geom.merge_geometry(stator.geom)
    geom.merge_geometry(rotor.geom)

    if di > 0:
        geom.add_line(Line(Node(0.0, 0.0), Node(di, 0.0)))
        geom.add_line(Line(Node(67.32, 0.0), Node(di + 67.32, 0.0)))

    start_1 = Node(2.0, 0.75)
    c1_1 = Node(X[0], X[1])  # it comes from X
    c2_1 = Node(X[2], X[3])  # it comes from X
    end_1 = Node(20.44, 0.75)
    geom.add_cubic_bezier(CubicBezier(start_1, c1_1, c2_1, end_1))

    start_1 = Node(24.44, 0.75)
    c1_1 = Node(X[4], X[5])  # it comes from X
    c2_1 = Node(X[6], X[7])  # it comes from X
    end_1 = Node(42.88, 0.75)
    geom.add_cubic_bezier(CubicBezier(start_1, c1_1, c2_1, end_1))

    start_1 = Node(46.88, 0.75)
    c1_1 = Node(X[8], X[9])  # it comes from X
    c2_1 = Node(X[10], X[11])  # it comes from X
    end_1 = Node(65.32, 0.75)
    geom.add_cubic_bezier(CubicBezier(start_1, c1_1, c2_1, end_1))

    geom.generate_intersections()
    # geom.export_svg(export_loc / f"geom-{int(di * 1000)}.svg")

    snapshot.add_geometry(geom)

    snapshot.assign_boundary_condition(0, 85, name="PB1")
    snapshot.assign_boundary_condition(67.33, 85, name="PB1")

    snapshot.assign_boundary_condition(0, 35, name="PB2")
    snapshot.assign_boundary_condition(67.33, 35, name="PB2")

    snapshot.assign_boundary_condition(0, 0.5, name="PB3")
    snapshot.assign_boundary_condition(67.33, 0.5, name="PB3")

    if di > 0:
        snapshot.assign_boundary_condition(di * 0.5, 0, name="PB4")
        snapshot.assign_boundary_condition(di * 0.5 + 67.33, 0, name="PB4")

    snapshot.assign_boundary_condition(di, -0.5, name="PB5")
    snapshot.assign_boundary_condition(di + 67.33, -0.45, name="PB5")

    snapshot.assign_boundary_condition(di, -6, name="PB6")
    snapshot.assign_boundary_condition(di + 67.33, -6, name="PB6")

    snapshot.assign_boundary_condition(di, -35, name="PB7")
    snapshot.assign_boundary_condition(di + 67, -35, name="PB7")

    snapshot.assign_boundary_condition(35, 100, name="a0")
    snapshot.assign_boundary_condition(di + 20, -50, name="a0")

    snapshot.add_postprocessing("integration", [(di + 10, -20), (di + 2, -6), (di + 65, -6), (di + 32, -5)], "Fx")

    snapshot.add_postprocessing("saveimage", basepath / "media" / "gif" / f"step-{int(di * 1000)}", None)

    return snapshot


if __name__ == "__main__":
    from random import uniform

    stator = ModelPiece("stator")
    stator.load_piece_from_dxf(basepath / "stator_stripped.dxf")
    stator.put(0, 0)
    rotor = ModelPiece("rotor")
    rotor.load_piece_from_dxf(basepath / "rotor.dxf")

    bounds = (
        (2, 12),
        (0.1, 1),
        (6, 20.44),
        (0.1, 1),
        (24.44, 35.5),
        (0.1, 1),
        (30, 42.88),
        (0.1, 1),
        (46.88, 60),
        (0.1, 1),
        (52, 65.32),
        (0.1, 1),
    )

    X = [uniform(lower, upper) for lower, upper in bounds]

    snapshot = generate_snapshot(stator, rotor, X, 4)

    snapshot.export()
    t0 = perf_counter()
    snapshot.execute()
    t1 = perf_counter()
    res = snapshot.retrive_results()
    print(res["Fx"])
    print("computation time:", t1 - t0)
