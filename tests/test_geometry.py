import os
from unittest import TestCase

from importlib_resources import files

from adze_modeler.geometry import Geometry

# from adze_modeler.gmsh import gmsh_writer
from adze_modeler.objects import CubicBezier, Line, Node, CircleArc


class TestGeometry(TestCase):
    def test_initialization(self):
        geo = Geometry()

        self.assertEqual([], geo.nodes)
        self.assertEqual([], geo.lines)
        self.assertEqual([], geo.cubic_beziers)

    def test_add_nodes_and_line(self):
        geo = Geometry()

        a = Node(1.0, 0.0)
        b = Node(0.5, 0.0)

        l = Line(a, b, id=1, label="test")

        geo.add_node(a)
        geo.add_node(b)
        self.assertEqual(a, geo.nodes[0])

        geo.add_line(l)
        self.assertEqual(l, geo.lines[0])

    def test_add_bezier(self):
        geo = Geometry()

        a = Node(1.0, 0.0)
        b = Node(0.5, 0.0)

        # l = Line(a, b, id=1, label="test")

        c1 = Node(0.6, 0.1)
        c2 = Node(0.7, 0.2)

        cb = CubicBezier(a, c1, c2, b, id=1, label="test")

        geo.add_node(a)
        geo.add_node(b)
        geo.add_node(c1)
        geo.add_node(c2)

        self.assertEqual(c1, geo.nodes[2])

        geo.add_cubic_bezier(cb)
        self.assertEqual(cb, geo.cubic_beziers[0])

    def test_intersections(self):
        path_rect = files("tests.svgtests").joinpath("rect.svg")
        path_overlap = files("tests.svgtests").joinpath("overlap.svg")

        g = Geometry()
        g.import_svg(str(path_rect))
        g.generate_intersections()

        self.assertEqual(len(g.nodes), 35)
        self.assertEqual(len(g.lines), 41)

        g = Geometry()
        g.import_svg(str(path_overlap))
        g.generate_intersections()

        self.assertEqual(len(g.nodes), 34)
        self.assertEqual(len(g.lines), 52)

    def test_append_node(self):

        geo = Geometry()
        a = Node(1.0, 0.0, id=1)
        b = Node(0.5, 0.0, id=2)
        c = Node(0.5000000000001, 0.0, id=3)

        geo.add_node(a)
        geo.add_node(b)
        self.assertEqual(2, len(geo.nodes))

        res = geo.append_node(c)
        self.assertEqual(b, res)
        self.assertEqual(2, len(geo.nodes))

    def test_merge_points(self):
        # after merging the nodes the number of them cannot increased

        geo = Geometry()

        a = Node(1.0, 0.0, id=1)
        b = Node(0.5, 0.0, id=2)
        c = Node(0.50000001, 0.0, id=3)
        d = Node(0.75, 0.75, id=6)

        l1 = Line(a, b, id=4, label="test1")
        l2 = Line(a, c, id=5, label="test2")
        l3 = Line(c, d, id=9, label="test3")

        geo.add_node(a)
        geo.add_node(b)
        geo.add_node(c)

        geo.add_line(l1)
        geo.add_line(l2)
        geo.add_line(l3)

        geo.merge_points()
        geo.merge_lines()

        self.assertEqual(3, len(geo.nodes))
        self.assertEqual(2, len(geo.lines))

        # calculate the total number of the nodes
        node_set = set()
        for i in geo.nodes:
            node_set.add(i.id)
        for l in geo.lines:
            node_set.add(l.start_pt.id)
            node_set.add(l.end_pt.id)

        self.assertEqual(len(node_set), 3)
        print(node_set)

    def test_find_surface(self):
        # in this test example a simple puzzle piece is defined by only simple and connected lines
        eml = files("tests.pygmsh_tests.test_cases").joinpath("test_lines.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())
        # set the tolerance to merge the given lines
        geo.epsilon = 1e-6
        geo.merge_points()
        geo.merge_lines()
        # there is only one described surface exists in the given geometry
        surfaces = geo.find_surfaces()

        self.assertEqual(len(surfaces), 1)

    @staticmethod
    def add_triangular_geometry():
        # creates a triangular shaped geometry which contains a line, a bezier line and a circle arc
        geo = Geometry()

        # points
        a = Node(1.0, 0.0)
        b = Node(0.0, 1.0)
        c = Node(-1.0, 0.0)

        origo = Node(0.0, 0.0)

        # control points for the bezier curve
        c1 = Node(-0.25, -0.9)
        c2 = Node(-0.75, -0.3)

        # adding the nodes
        geo.add_node(a)
        geo.add_node(b)
        geo.add_node(c)

        line = Line(a,b)
        bezier = CubicBezier(start_pt=b,  control1=c1, control2=c2, end_pt=c)
        circle_arc = CircleArc(c,origo,a)

        print(bezier)

        geo.add_line(line)
        geo.add_arc(circle_arc)
        geo.add_cubic_bezier(bezier)
        return geo

    def test_find_surface_loop(self):
        # in the case of lines, bezier curves and circle arcs
        geo = self.add_triangular_geometry()

        # adding lines, beziers and circle arcs cannot create new nodes
        self.assertEqual(len(geo.nodes), 3)
        self.assertEqual(len(geo.circle_arcs),1)
        self.assertEqual(len(geo.cubic_beziers),1)
        self.assertEqual(len(geo.lines),1)
        geo.export_svg()
        return
