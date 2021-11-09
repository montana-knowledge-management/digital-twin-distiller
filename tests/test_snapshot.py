import unittest
from pathlib import Path

from digital_twin_distiller.boundaries import DirichletBoundaryCondition, NeumannBoundaryCondition
from digital_twin_distiller.geometry import Geometry
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import Agros2DMetadata, FemmMetadata
from digital_twin_distiller.objects import Line, Node
from digital_twin_distiller.platforms.agros2d import Agros2D
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot


class MockFileHandle:
    def __init__(self):
        self.content = ""

    def close(self):
        pass

    def write(self, str_):
        self.content += str_

    def clear(self):
        self.content = ""

    def get_line(self, n: int):
        return self.content.split("\n")[n]


class TestSnapshotAgros2D(unittest.TestCase):
    def get_geometry(self):
        g = Geometry()
        g.add_line(Line(Node(-1, 0), Node(1, 0)))
        return g

    def get_metadata(self):
        agros_metadata = Agros2DMetadata()
        agros_metadata.file_script_name = "agros_solver_script"
        agros_metadata.file_metrics_name = "agros_solution.csv"
        agros_metadata.problem_type = "magnetic"
        agros_metadata.coordinate_type = "axisymmetric"
        agros_metadata.analysis_type = "steadystate"
        agros_metadata.unit = 1e-3
        agros_metadata.nb_refinements = 2
        return agros_metadata

    def get_platform(self):
        return Agros2D(self.get_metadata())

    def get_snapshot(self):
        return Snapshot(self.get_platform())

    def test_snapshot_setup(self):
        metadata = self.get_metadata()
        platform = self.get_platform()
        snapshot = Snapshot(platform)

        snapshot.set_platform(platform)
        self.assertTrue(
            snapshot.platform.metadata.compatible_platform,
            metadata.compatible_platform,
        )

    def test_set_add_boundary_condition(self):
        s = self.get_snapshot()
        s.add_boundary_condition(DirichletBoundaryCondition("eper", "magnetic", magnetic_potential=3))
        self.assertTrue("eper" in s.boundaries)

    def test_set_add_neumann_boundary_condition(self):
        s = self.get_snapshot()
        s.add_boundary_condition(NeumannBoundaryCondition("malna", "magnetic", surface_current=1.5))
        self.assertTrue("malna" in s.boundaries)

    def test_assign_boundary_condition(self):
        s = self.get_snapshot()
        s.add_boundary_condition(DirichletBoundaryCondition("eper", "magnetic", magnetic_potential=3))
        self.assertRaises(ValueError, s.assign_boundary_condition, x=0, y=0, name="falsename")

    def test_add_material(self):
        s = self.get_snapshot()
        s.add_material(Material("iron"))
        self.assertTrue(s.materials)

    def test_add_material(self):
        s = self.get_snapshot()
        s.add_material(Material("iron"))
        s.assign_material(4, 42, "iron")
        self.assertEqual(s.materials["iron"].assigned[0], (4, 42))

    def test_add_geometry(self):
        g = self.get_geometry()
        s = self.get_snapshot()
        s.add_geometry(g)
        self.assertTrue(len(s.lines) == 1)
        self.assertTrue(len(s.nodes) == 2)

    def test_export_results(self):
        s = self.get_snapshot()
        s.add_postprocessing('mesh_info', None, None)
        s.add_postprocessing('point_value', [1.0, 1.0], "Bx")
        s.add_postprocessing('integration', [(1.0, 2.0), (3.0, 3.0)], "Energy")
        f = MockFileHandle()
        s.export(f)

        # mesh data
        self.assertIn(r'info = magnetic.solution_mesh_info()', f.content)
        self.assertIn(r'f.write("{}, {}\n".format("dofs", info["dofs"]))', f.content)
        self.assertIn(r'f.write("{}, {}\n".format("nodes", info["nodes"]))', f.content)
        self.assertIn(r'f.write("{}, {}\n".format("elements", info["elements"])', f.content)

        # point value
        self.assertIn(r'point = magnetic.local_values(0.001, 0.001)["Brx"]', f.content)
        self.assertIn(r'f.write("{}, 0.001, 0.001, {}\n".format("Bx", point))', f.content)

        # integral value
        self.assertIn(r"magnetic.volume_integrals([(1.0, 2.0), (3.0, 3.0)])", f.content)
        self.assertIn(r'f.write("Energy, {}\n".format(val))', f.content)

    # TODO: ??
    def test_export(self):
        # compares the created model with a reference script
        s = self.get_snapshot()
        s.add_material(Material("air"))
        s.assign_material(0, 0, "air")
        s.add_geometry(self.get_geometry())
        s.add_boundary_condition(DirichletBoundaryCondition("d0", field_type="magnetic", magnetic_potential=30))
        s.add_boundary_condition(NeumannBoundaryCondition("n0", field_type="magnetic", surface_current=0))
        f = MockFileHandle()
        # s.export()
        s.export(f)
        with open(Path(__file__).parent / "solver_script_references" / "agros2d_default_script.py") as f_ref:
            for i, line in enumerate(f_ref.readlines()):
                if "remove(" in line:
                    continue

                if "openfile" in line:
                    continue

                if "saveas" in line:
                    continue

                self.assertEqual(f.get_line(i), line.rstrip())


class TestSnapshotFemm(unittest.TestCase):
    def get_geometry(self):
        g = Geometry()
        g.add_line(Line(Node(-1, 0), Node(1, 0)))
        return g

    def get_metadata(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "axisymmetric"
        femm_metadata.file_script_name = "femm_solver_script"
        femm_metadata.file_metrics_name = "femm_solution.csv"
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False
        return femm_metadata

    def get_platform(self):
        return Femm(self.get_metadata())

    def get_snapshot(self):
        return Snapshot(self.get_platform())

    def test_snapshot_setup(self):
        metadata = self.get_metadata()
        platform = self.get_platform()
        snapshot = Snapshot(platform)

        snapshot.set_platform(platform)
        self.assertTrue(
            snapshot.platform.metadata.compatible_platform,
            metadata.compatible_platform,
        )

    def test_set_add_boundary_condition(self):
        s = self.get_snapshot()
        s.add_boundary_condition(DirichletBoundaryCondition("eper", "magnetic", magnetic_potential=3))
        self.assertTrue("eper" in s.boundaries)

    def test_assign_boundary_condition(self):
        s = self.get_snapshot()
        s.add_boundary_condition(DirichletBoundaryCondition("eper", "magnetic", magnetic_potential=3))
        self.assertRaises(ValueError, s.assign_boundary_condition, x=0, y=0, name="falsename")

    def test_add_material(self):
        s = self.get_snapshot()
        s.add_material(Material("iron"))
        self.assertTrue(s.materials)

    def test_add_material(self):
        s = self.get_snapshot()
        s.add_material(Material("iron"))
        s.assign_material(4, 42, "iron")
        self.assertEqual(s.materials["iron"].assigned[0], (4, 42))

    def test_add_geometry(self):
        g = self.get_geometry()
        s = self.get_snapshot()
        s.add_geometry(g)
        self.assertTrue(len(s.lines) == 1)
        self.assertTrue(len(s.nodes) == 2)

    def test_export_results(self):
        s = self.get_snapshot()
        s.add_postprocessing('mesh_info', None, None)
        s.add_postprocessing('point_value', [1.0, 1.0], "Bx")
        s.add_postprocessing('integration', [(1.0, 2.0), (3.0, 3.0)], "Energy")

        f = MockFileHandle()
        s.export(f)

        # mesh info
        self.assertIn(r'write(file_out, "nodes, ", mo_numnodes(), "\n")', f.content)
        self.assertIn(r'write(file_out, "elements, ", mo_numelements(), "\n")', f.content)

        # point values
        self.assertIn(r'A, B1, B2, Sig, E, H1, H2, Je, Js, Mu1, Mu2, Pe, Ph = mo_getpointvalues(1.0, 1.0)', f.content)
        self.assertIn(r'write(file_out, "Bx, 1.0, 1.0, ", B1, "\n")', f.content)

        # energy
        self.assertIn(r'mo_selectblock(1.0, 2.0)', f.content)
        self.assertIn(r'Energy = mo_blockintegral(2)', f.content)

    def test_export(self):
        s = self.get_snapshot()
        s.add_geometry(self.get_geometry())
        s.add_material(Material("air"))
        s.assign_material(0, 0, "air")
        s.add_boundary_condition(DirichletBoundaryCondition("d0", "magnetic", magnetic_potential=30))
        f = MockFileHandle()
        # s.export()
        s.export(f)
        with open(Path(__file__).parent / "solver_script_references" / "femm_default_script.lua") as f_ref:
            for i, line in enumerate(f_ref.readlines()):
                if "remove(" in line:
                    continue

                if "openfile" in line:
                    continue

                if "saveas" in line:
                    continue

                self.assertEqual(f.get_line(i), line.rstrip())


if __name__ == "__main__":
    t = TestSnapshotAgros2D()
    t.test_export()
