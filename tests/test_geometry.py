import os
from unittest import TestCase

from importlib_resources import files

from adze_modeler.geometry import Geometry

# from adze_modeler.gmsh import gmsh_writer
from adze_modeler.objects import CubicBezier, Line, Node


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


# class TestMeshing(TestCase):
#     def test_mesh_the_triangle(self):
#         path = files("examples.triangle").joinpath("triangle.svg")
#         # print(path)
#         geo = Geometry()
#         geo.import_svg(path.as_posix())
#
#         node = Node(555, -555)
#         geo.add_node(node)
#
#         res = geo.__str__()
#         self.assertIn(node.__str__(), res)
#         print(res)
#         gmsh_writer(geo.nodes, geo.lines, geo.circle_arcs, geo.cubic_beziers)
#
#         try:
#             os.remove("test.geo_unrolled")
#         except FileNotFoundError:
#             try:
#                 # running tox
#                 os.remove("tests\test.geo_unrolled")
#             except FileNotFoundError:
#                 self.assertTrue(False)
#
#     def test_mesh_the_owl(self):
#         path = files("examples.owl").joinpath("owl-shape.svg")
#         geo = Geometry()
#         geo.import_svg(path.as_posix())
#         gmsh_writer(geo.nodes, geo.lines, geo.circle_arcs, geo.cubic_beziers)
#
#         try:
#             os.remove("test.geo_unrolled")
#         except FileNotFoundError:
#             try:
#                 # running tox
#                 os.remove("tests\test.geo_unrolled")
#             except FileNotFoundError:
#                 self.assertTrue(False)
