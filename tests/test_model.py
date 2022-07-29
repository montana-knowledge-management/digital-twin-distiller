import unittest
from pathlib import Path
from shutil import rmtree

from digital_twin_distiller import DirichletBoundaryCondition, Line, Material, Node
from digital_twin_distiller.metadata import Agros2DMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.platforms.agros2d import Agros2D
from digital_twin_distiller.snapshot import Snapshot


def mock_export(*args, **kwargs):
    return True


def mock_execute(*args, **kwargs):
    return True


def mock_results():
    return {"a": 3, "b": 4}


class MockModel(BaseModel):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def setup_solver(self):
        super().setup_solver()
        # Agros2D
        agros_metadata = Agros2DMetadata()
        agros_metadata.file_script_name = self.file_solver_script
        agros_metadata.file_metrics_name = self.file_solution
        agros_metadata.problem_type = "magnetic"
        agros_metadata.coordinate_type = "axisymmetric"
        agros_metadata.analysis_type = "steadystate"
        agros_metadata.unit = 1e-3
        agros_metadata.nb_refinements = 0
        agros_metadata.adaptivity = "hp-adaptivity"
        agros_metadata.adaptivity_tol = 1
        platform_agros = Agros2D(agros_metadata)
        self.snapshot = Snapshot(platform_agros)

        self.snapshot.export = mock_export
        self.snapshot.execute = mock_execute
        self.snapshot.retrive_results = mock_results

    def add_postprocessing(self):
        self.snapshot.add_postprocessing("mesh_info", None, None)

    def define_materials(self):
        m = Material("testmaterial")
        self.snapshot.add_material(m)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("testboundary", "magnetic")

        self.snapshot.add_boundary_condition(a0)

    def build_geometry(self):
        self.snapshot.add_geometry(self.geom)


class TestModel(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dir_current = Path(__file__).parent
        self.dir_resources = self.dir_current / "resources"
        self.dir_snapshots = self.dir_current / "snapshots"
        self.dir_media = self.dir_current / "media"
        self.dir_data = self.dir_current / "data"
        # self.dir_export = self.dir_snapshots / self.name
        self.dir_export = self.dir_snapshots

    def test_model_paths(self):
        m = MockModel(exportname="test_name")
        self.assertEqual(self.dir_current, m.dir_current)
        self.assertEqual(self.dir_resources, m.dir_resources)
        self.assertEqual(self.dir_snapshots, m.dir_snapshots)
        self.assertEqual(self.dir_media, m.dir_media)
        self.assertEqual(self.dir_data, m.dir_data)
        self.assertEqual(self.dir_export, m.dir_export.parent)

        # create paths
        m._init_directories()
        self.assertTrue(self.dir_current.exists())
        self.assertTrue(self.dir_resources.exists())
        self.assertTrue(self.dir_snapshots.exists())
        self.assertTrue(self.dir_media.exists())
        self.assertTrue(self.dir_data.exists())
        self.assertTrue(self.dir_export.exists())

        rmtree(m.dir_snapshots)

        # teardown
        self.dir_resources.rmdir()
        self.dir_media.rmdir()
        self.dir_data.rmdir()

    def test_add_line(self):
        m = MockModel(exportname="test_name")

        self.assertFalse(m.geom.nodes)
        self.assertFalse(m.geom.lines)
        self.assertFalse(m.geom.circle_arcs)
        self.assertFalse(m.geom.cubic_beziers)

        n0 = Node(0, 0)
        n1 = Node(5.3, -3.2)
        l1 = Line(n0, n1)
        m.add_line(0, 0, 5.3, -3.2)

        self.assertEqual(len(m.geom.nodes), 2)
        self.assertEqual(len(m.geom.lines), 1)
        self.assertFalse(m.geom.circle_arcs)
        self.assertFalse(m.geom.cubic_beziers)

        self.assertEqual(m.geom.lines[0], l1)

    def test_add_circle_arc(self):
        m = MockModel(exportname="test_name")

        self.assertFalse(m.geom.nodes)
        self.assertFalse(m.geom.lines)
        self.assertFalse(m.geom.circle_arcs)
        self.assertFalse(m.geom.cubic_beziers)

        m.add_circle_arc(-1, 0, 0, 0, 1, 0)

        self.assertEqual(len(m.geom.nodes), 2)
        self.assertFalse(m.geom.lines)
        self.assertEqual(len(m.geom.circle_arcs), 1)
        self.assertFalse(m.geom.cubic_beziers)

        self.assertTrue(Node(-1, 0) in m.geom.nodes)
        self.assertTrue(Node(1, 0) in m.geom.nodes)
        self.assertFalse(Node(0, 0) in m.geom.nodes)
        self.assertFalse(Node(0, -1) in m.geom.nodes)

    def test_assign_material(self):
        m = MockModel(exportname="test_name")

        self.assertFalse(m.label_queue)

        m.assign_material(0, 1, "testmaterial")

        self.assertEqual(len(m.label_queue), 1)
        self.assertEqual(m.label_queue[0], (0, 1, "testmaterial"))

    def test_assign_boundary(self):
        m = MockModel(exportname="test_name")

        self.assertFalse(m.label_queue)

        m.assign_boundary(0, 1, "testboundary")

        self.assertEqual(len(m.boundary_queue), 1)
        self.assertEqual(m.boundary_queue[0], (0, 1, "testboundary"))

    def test_assign_boundary_arc(self):
        m = MockModel(exportname="test_name")

        self.assertFalse(m.label_queue)

        m.assign_boundary_arc(0, 1, "testboundary")

        self.assertEqual(len(m.boundary_arc_queue), 1)
        self.assertEqual(m.boundary_arc_queue[0], (0, 1, "testboundary"))

    def test_call(self):

        m = MockModel(exportname="test_name")
        m.add_line(1, 1, 1, -1)
        m.add_circle_arc(-1, 0, 0, 0, 1, 0)
        m.assign_material(0, 1, "testmaterial")
        m.assign_boundary(1, 0, "testboundary")
        m.assign_boundary_arc(0, 1, "testboundary")

        res = m(cleanup=False, devmode=False)
        self.assertTrue(res["a"], 5)
        self.assertTrue(res["b"], 6)

        src = Path(__file__).parent / "snapshots/test_name"
        src.mkdir(parents=True, exist_ok=True)
        with open(src / "P_test_name.py", "w") as f:
            f.write("testtext")
        res = m(cleanup=True, devmode=False)
        self.assertTrue(res["a"], 5)
        self.assertTrue(res["b"], 6)

        # this should raise an exception therefore
        # the return value should be a None
        def f():
            raise

        m.build_geometry = f
        res = m(cleanup=False, devmode=False)
        self.assertTrue(res is None)
