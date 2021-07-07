from pathlib import Path
from numpy import linspace

from adze_modeler.geometry import Geometry
from adze_modeler.metadata import FemmMetadata
from adze_modeler.objects import Node, Line, CubicBezier
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.modelpiece import ModelPiece

basepath = Path(__file__).parent
exportpath = basepath / "snapshots"
exportpath.mkdir(exist_ok=True)

femm_metadata = FemmMetadata()
femm_metadata.problem_type = "magnetic"
femm_metadata.coordinate_type = "planar"
femm_metadata.file_script_name = basepath / "femm_solver_script"
femm_metadata.file_metrics_name = basepath / "femm_solution.csv"
femm_metadata.unit = "millimeters"
femm_metadata.smartmesh = False
platform = Femm(femm_metadata)

stator = ModelPiece("stator")
stator.load_piece_from_dxf(basepath / "stator_stripped.dxf")
stator.put(0, 0)
rotor = ModelPiece("rotor")
rotor.load_piece_from_dxf(basepath / "rotor.dxf")
for i, di in enumerate(linspace(0, 50, 3)):
    rotor.put(di, -50.75)
    geom = Geometry()
    geom.merge_geometry(stator.geom)
    geom.merge_geometry(rotor.geom)
    if i > 0:
        geom.add_line(Line(Node(0.0, 0.0), Node(di, 0.0)))
        geom.add_line(Line(Node(67.32, 0.0), Node(di + 67.32, 0.0)))

    start_1 = Node(2.0, 0.75)
    c1_1 = Node(2.5, 0) # it comes from X
    c2_1 = Node(3.5, 0) # it comes from X
    end_1 = Node(20.44, 0.75)
    geom.add_cubic_bezier(CubicBezier(start_1, c1_1, c2_1, end_1))

    start_1 = Node(24.44, 0.75)
    c1_1 = Node(30, 0)  # it comes from X
    c2_1 = Node(35, 0)  # it comes from X
    end_1 = Node(42.88, 0.75)
    geom.add_cubic_bezier(CubicBezier(start_1, c1_1, c2_1, end_1))

    start_1 = Node(46.88, 0.75)
    c1_1 = Node(50, 0)  # it comes from X
    c2_1 = Node(60, 0)  # it comes from X
    end_1 = Node(65.32, 0.75)
    geom.add_cubic_bezier(CubicBezier(start_1, c1_1, c2_1, end_1))

    geom.generate_intersections()
    geom.export_svg(exportpath / f"geom-{i}.svg")

snapshot = Snapshot(platform)
snapshot.add_geometry(geom)