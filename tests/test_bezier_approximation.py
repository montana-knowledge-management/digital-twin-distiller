import unittest
import uuid
from adze_modeler.geometry import Geometry
from adze_modeler.objects import CubicBezier
from adze_modeler.objects import Line
from adze_modeler.objects import Node

from importlib_resources import files


class TestBezierApprox(unittest.TestCase):
    def test_casteljeu_step(self):
        # Test Nodes
        a = Node(x=0.0, y=0.0, id=1)
        b = Node(x=1.0, y=-0.2, id=2)
        c1 = Node(x=0.2, y=0.9, id=3)
        c2 = Node(x=0.8, y=0.7, id=4)

        bez = CubicBezier(a, c1, c2, b)

        r, l = Geometry.casteljau(bez)

        self.assertEqual(r.end_pt.x, 1.0)
        self.assertEqual(r.end_pt.y, -0.2)

        self.assertEqual(l.control1.x, 0.1)
        self.assertEqual(l.control1.y, 0.45)

    # def test_1(self):
    #     #eml = files("tests.svgtests").joinpath("test_bez_approx.svg")
    #     eml = files("examples.owl").joinpath("owl-shape.svg")
    #
    #     geo = Geometry()
    #     geo.import_svg(eml.as_posix())
    #
    #     for bezier in geo.cubic_beziers:
    #         r, l = geo.casteljau(bezier)
    #         geo.add_line(Line(r.start_pt, r.control1))
    #         geo.add_line(Line(r.control1, r.control2))
    #         geo.add_line(Line(r.control2, r.end_pt))
    #
    #         geo.add_line(Line(l.start_pt, l.control1))
    #         geo.add_line(Line(l.control1, l.control2))
    #         geo.add_line(Line(l.control2, l.end_pt))
    #
    #     print(geo.lines)
    #     print(len(geo.lines))
    #
    #     geo.export_svg('test.svg')
