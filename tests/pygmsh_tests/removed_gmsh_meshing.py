from os import remove
from unittest import TestCase

from importlib_resources import files
from meshio._helpers import read

from adze_modeler.geometry import Geometry
from adze_modeler.gmsh import GMSHModel
from adze_modeler.objects import CircleArc, Line, Node

# plotting out the mesh
import pyvista as pv

# msh = pv.read(file_name + '.vtk')
# msh.plot(show_edges=True)


class TestGMSHWriter(TestCase):
    def test_only_line_surface(self):
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

        # create a gmsh mesh from the given geometry
        gmsh = GMSHModel(geo)
        gmsh.gmsh_writer("test1")

        msh_data = read("test1.vtk")

        # tests that the code generates a valid mesh
        self.assertGreaterEqual(len(msh_data.cells), 3)

        # check the surface, the surface should contain 8 edges
        self.assertEqual(len(surfaces[0]), 8)
        self.assertEqual(round(surfaces[0][0].start_pt.x, 1), 103.4)
        # remove the geo and msh files
        remove("test1.geo_unrolled")
        remove("test1.vtk")
#
#     def test_bezier_line_surface(self):
#         # import the surface
#         eml = files("tests.pygmsh_tests.test_cases").joinpath("test_bezier.svg")
#         geo = Geometry()
#         geo.import_svg(eml.as_posix())
#         # set the tolerance to merge the given lines
#         geo.epsilon = 1e-6
#         geo.merge_points()
#
#         # there is only one described surface exists in the given geometry
#         surfaces = geo.find_surfaces()
#         geo.plot_connection_graph()
#         # create a gmsh mesh from the given geometry
#         gmsh = GMSHModel(geo)
#         gmsh.gmsh_writer("test2")
#
#         msh_data = read("test2.vtk")
#
#         # tests that the code generates a valid mesh
#         self.assertGreaterEqual(len(msh_data.cells), 3)
#
#         # check the surface, the surface should contain only 7 edges
#         self.assertEqual(len(surfaces[0]), 7)
#         self.assertEqual(round(surfaces[0][0].start_pt.x, 1), 100.1)
#
#         # remove the geo and msh files
#         remove("test2.geo_unrolled")
#         remove("test2.vtk")
#
#     def test_circle_defined_surface(self):
#         # define the geometry by hand, a simple arc
#         geo = Geometry()
#
#         a = Node(x=0.0, y=0.0, id=1)
#         b = Node(x=10.0, y=0.0, id=2)
#         c = Node(x=0.0, y=10.0, id=3)
#
#         geo.add_line(Line(a, b))
#         geo.add_line(Line(a, c))
#         geo.add_arc(CircleArc(c, a, b))
#
#         # there is only one described surface exists in the given geometry
#         surfaces = geo.find_surfaces()
#
#         gmsh = GMSHModel(geo)
#         gmsh.gmsh_writer("test3")
#
#         # check the surface, the surface should contain only 3 edges
#         self.assertEqual(len(surfaces[0]), 3)
#         self.assertEqual(round(surfaces[0][0].start_pt.x, 1), 0.0)
#
#         # remove the geo and msh files
#         remove("test3.geo_unrolled")
#         remove("test3.vtk")
