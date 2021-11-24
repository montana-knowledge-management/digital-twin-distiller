from io import StringIO
from unittest import TestCase

from digital_twin_distiller import (
    DirichletBoundaryCondition,
    Geometry,
    Line,
    Material,
    NgElectrostaticMetadata,
    NgElectrostatics,
    Node,
)
from digital_twin_distiller.snapshot import Snapshot

"""
This file is for testing the NgSolve base class and its ability to
detect surfaces
"""


class TestSurfaceDetection(TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        ng_metadata = NgElectrostaticMetadata()
        ng_metadata.file_script_name = "testscriptname"
        platform = NgElectrostatics(ng_metadata)

        cls.snapshot = Snapshot(platform)

        mr2 = Material("r2")
        mr3 = Material("r3")
        mr4 = Material("r4")
        mr5 = Material("r5")
        mr6 = Material("r6")

        cls.snapshot.add_material(mr2)
        cls.snapshot.add_material(mr3)
        cls.snapshot.add_material(mr4)
        cls.snapshot.add_material(mr5)
        cls.snapshot.add_material(mr6)

        a0 = DirichletBoundaryCondition("a0", field_type="electrostatic", fixed_voltage=0.0)
        cls.snapshot.add_boundary_condition(a0)

        points = [(0, 0)]
        cls.snapshot.add_postprocessing("integration", points, "Energy")

        cls.p1 = Node(4.3, 3.68, label="p1")
        cls.p2 = Node(6.38, 5.57, label="p2")
        cls.p3 = Node(10.39, 6.69, label="p3")
        cls.p4 = Node(14.03, 6.09, label="p4")
        cls.p5 = Node(13.75, 2.34, label="p5")
        cls.p6 = Node(11.33, 0.23, label="p6")
        cls.p7 = Node(6.57, 1.22, label="p7")
        cls.p8 = Node(6.41, 4.16, label="p8")
        cls.p9 = Node(8.04, 3.3, label="p9")
        cls.p10 = Node(8.79, 4.35, label="p10")
        cls.p11 = Node(12.62, 3.34, label="p11")

        g = Geometry()

        g.add_line(Line(cls.p1, cls.p2))
        g.add_line(Line(cls.p2, cls.p3))
        g.add_line(Line(cls.p3, cls.p4))
        g.add_line(Line(cls.p4, cls.p5))
        g.add_line(Line(cls.p5, cls.p6))
        g.add_line(Line(cls.p6, cls.p7))
        g.add_line(Line(cls.p7, cls.p1))
        g.add_line(Line(cls.p7, cls.p11))
        g.add_line(Line(cls.p11, cls.p5))
        g.add_line(Line(cls.p11, cls.p3))
        g.add_line(Line(cls.p10, cls.p2))
        g.add_line(Line(cls.p10, cls.p9))
        g.add_line(Line(cls.p9, cls.p8))
        g.add_line(Line(cls.p8, cls.p1))
        g.add_line(Line(cls.p10, cls.p3))

        cls.snapshot.add_geometry(g)

        cls.snapshot.assign_material(8, 5, "r2")
        cls.snapshot.assign_material(6.5, 5.2, "r3")
        cls.snapshot.assign_material(7.5, 3, "r4")
        cls.snapshot.assign_material(10, 1.3, "r5")
        cls.snapshot.assign_material(12, 6, "r6")

    @classmethod
    def tearDownClass(cls) -> None:
        ...

    def test_surfaces(self):
        output = StringIO()
        self.snapshot.export(output)
        # get the surfaces
        surfaces = self.snapshot.platform.compose_geometry()

        # convert the surface into sets
        surfaces = [set(si.nodes) for si in surfaces]

        r2 = {self.p3, self.p2, self.p10}
        r3 = {self.p2, self.p1, self.p8, self.p9, self.p10}
        r4 = {self.p1, self.p7, self.p11, self.p3, self.p10, self.p9, self.p8}
        r5 = {self.p7, self.p6, self.p5, self.p11}
        r6 = {self.p5, self.p11, self.p3, self.p4}

        assert len(surfaces) == 5
        assert any(r2 == si for si in surfaces)
        assert any(r3 == si for si in surfaces)
        assert any(r4 == si for si in surfaces)
        assert any(r5 == si for si in surfaces)
        assert any(r6 == si for si in surfaces)

    def test_render(self):
        # Its enough to pass if the method does not throw any exeption.
        output = StringIO()
        self.snapshot.export(output)
        self.snapshot.platform.compose_geometry()
        self.snapshot.platform.render_geo()

    def test_execution(self):
        output = StringIO()
        self.snapshot.export(output)
        self.snapshot.platform.execute(runner="python")
