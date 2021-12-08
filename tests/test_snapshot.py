import math
import unittest
from pathlib import Path

from digital_twin_distiller.boundaries import (
    AntiPeriodicAirGap,
    AntiPeriodicBoundaryCondition,
    DirichletBoundaryCondition,
    NeumannBoundaryCondition,
    PeriodicAirGap,
    PeriodicBoundaryCondition,
)
from digital_twin_distiller.femm_wrapper import femm_current_flow, femm_electrostatic, femm_heat_flow, femm_magnetic
from digital_twin_distiller.geometry import Geometry
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import Agros2DMetadata, FemmMetadata
from digital_twin_distiller.objects import CircleArc, Line, Node
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
        s.add_postprocessing("mesh_info", None, None)
        s.add_postprocessing("point_value", [1.0, 1.0], "Bx")
        s.add_postprocessing("integration", [(1.0, 2.0), (3.0, 3.0)], "Energy")
        f = MockFileHandle()
        s.export(f)

        # mesh data
        self.assertIn(r"info = magnetic.solution_mesh_info()", f.content)
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

        self.assertRaises(Exception, s.execute(cleanup=True))


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

    def test_platform_init(self):
        metadata = self.get_metadata()
        platform = Femm(metadata)
        self.assertEqual(platform.writer.field, femm_magnetic)

        metadata = self.get_metadata()
        metadata.problem_type = "electrostatic"
        platform = Femm(metadata)
        self.assertEqual(platform.writer.field, femm_electrostatic)

        metadata = self.get_metadata()
        metadata.problem_type = "heat"
        platform = Femm(metadata)
        self.assertEqual(platform.writer.field, femm_heat_flow)

        metadata = self.get_metadata()
        metadata.problem_type = "current"
        platform = Femm(metadata)
        self.assertEqual(platform.writer.field, femm_current_flow)

        metadata = self.get_metadata()
        metadata.problem_type = "falsefield"
        with self.assertRaises(ValueError):
            Femm(metadata)

    def test_set_add_boundary_condition(self):
        s = self.get_snapshot()
        s.add_boundary_condition(DirichletBoundaryCondition("eper", "magnetic", magnetic_potential=3))
        self.assertTrue("eper" in s.boundaries)

    def test_assign_boundary_condition(self):
        s = self.get_snapshot()
        # dirichlet
        s.add_boundary_condition(DirichletBoundaryCondition("eper", "magnetic", magnetic_potential=3))
        # neumann
        s.add_boundary_condition(NeumannBoundaryCondition("cekla", "magnetic", surface_current=12.4))
        # anti-periodic
        s.add_boundary_condition(AntiPeriodicBoundaryCondition("retek", "magnetic"))
        # periodic
        s.add_boundary_condition(AntiPeriodicBoundaryCondition("mogyoro", "magnetic"))

        # false boundary condition allocated
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
        s.add_postprocessing("mesh_info", None, None)
        s.add_postprocessing("point_value", [1.0, 1.0], "Bx")
        s.add_postprocessing("integration", [(1.0, 2.0), (3.0, 3.0)], "Energy")

        s.add_postprocessing("saveimage", None, None)

        f = MockFileHandle()
        s.export(f)

        # mesh info
        self.assertIn(r'write(file_out, "nodes, ", mo_numnodes(), "\n")', f.content)
        self.assertIn(r'write(file_out, "elements, ", mo_numelements(), "\n")', f.content)

        # point values
        self.assertIn(r"A, B1, B2, Sig, E, H1, H2, Je, Js, Mu1, Mu2, Pe, Ph = mo_getpointvalues(1.0, 1.0)", f.content)
        self.assertIn(r'write(file_out, "Bx, 1.0, 1.0, ", B1, "\n")', f.content)

        # energy
        self.assertIn(r"mo_selectblock(1.0, 2.0)", f.content)
        self.assertIn(r"Energy = mo_blockintegral(2)", f.content)

        # saveimage
        self.assertIn(r"mo_refreshview()", f.content)
        self.assertIn(r'mo_save_bitmap("', f.content)

    def test_export(self):
        s = self.get_snapshot()
        s.add_geometry(self.get_geometry())
        s.add_material(Material("air"))
        s.assign_material(0, 0, "air")
        s.add_boundary_condition(DirichletBoundaryCondition("d0", "magnetic", magnetic_potential=30))
        s.add_boundary_condition(NeumannBoundaryCondition("cekla", "magnetic", surface_current=12.4))
        s.add_boundary_condition(AntiPeriodicBoundaryCondition("retek", "magnetic"))
        s.add_boundary_condition(PeriodicBoundaryCondition("mogyoro", "magnetic"))
        s.add_boundary_condition(AntiPeriodicAirGap("g", "magnetic"))
        s.add_boundary_condition(PeriodicAirGap("f0", "magnetic"))

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

        # running should cause an exception
        self.assertRaises(Exception, s.execute(cleanup=True))

    def test_export_geometry_entities(self):
        s = self.get_snapshot()
        s.add_boundary_condition(DirichletBoundaryCondition("eper", "magnetic", magnetic_potential=3))
        g = Geometry()
        f = MockFileHandle()
        n0 = Node(0, 0)
        n1 = Node(1, 0)
        n2 = Node(1, 1)
        n3 = Node(0, 1)

        l1 = Line(n0, n1)
        l2 = Line(n1, n2)
        l3 = Line(n2, n3)
        l4 = Line(n3, n0)

        c1 = CircleArc(Node(0, 0.5), Node(0.5, 0.5), Node(1, 0.5))

        g.add_node(n0)
        g.add_node(n1)
        g.add_node(n2)
        g.add_node(n3)

        g.add_line(l1)
        g.add_line(l2)
        g.add_line(l3)
        g.add_line(l4)

        g.add_arc(c1)

        s.add_geometry(g)
        s.assign_boundary_condition(0.5, 0, name="eper")
        s.export(f)

        self.assertIn(r"mi_addnode(0, 0)", f.content)
        self.assertIn(r"mi_addnode(1, 0)", f.content)
        self.assertIn(r"mi_addnode(1, 1)", f.content)
        self.assertIn(r"mi_addnode(0, 1)", f.content)
        self.assertIn(r"mi_addnode(0, 0.5)", f.content)
        self.assertIn(r"mi_addnode(1, 0.5)", f.content)
        self.assertIn(r"mi_addsegment(0, 0, 1, 0)", f.content)
        self.assertIn(r"mi_selectsegment(0.5, 0.0)", f.content)
        self.assertIn(r'mi_setsegmentprop("eper", None, 1, 0, 0, "<None>")', f.content)
        self.assertIn(r"mi_clearselected()", f.content)
        self.assertIn(r"mi_addsegment(1, 0, 1, 1)", f.content)
        self.assertIn(r"mi_addsegment(1, 1, 0, 1)", f.content)
        self.assertIn(r"mi_addsegment(0, 1, 0, 0)", f.content)
        self.assertIn(r"mi_addarc(0, 0.5, 1, 0.5, 180.0, 20)", f.content)

    def test_retrive_results(self):
        s = self.get_snapshot()
        result_file = Path(s.platform.metadata.file_metrics_name)

        with open(result_file, "w", encoding="utf-8") as f:
            print("dofs, 10168", file=f)
            print("nodes, 2360", file=f)
            print("elements, 1125", file=f)
            print("Bx, 1, 2, 33", file=f)

        results = s.retrive_results()
        result_file.unlink()

        self.assertEqual(10168, results["dofs"])
        self.assertEqual(2360, results["nodes"])
        self.assertEqual(1125, results["elements"])

        self.assertTrue(len(results["Bx"]))
        self.assertListEqual([1.0, 2.0, 33.0], list(results["Bx"][0]))
