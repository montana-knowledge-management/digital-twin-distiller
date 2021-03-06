from os import remove
from unittest import TestCase

from importlib_resources import files
from meshio._helpers import read

from digital_twin_distiller.geometry import Geometry
from digital_twin_distiller.gmsh import GMSHModel
from digital_twin_distiller.objects import CircleArc, Line, Node


class TestGMSHWriter(TestCase):
    def test_only_line_surface(self):
        # in this test example a simple puzzle piece is defined by only simple and connected lines
        eml = files("tests.pygmsh_tests.test_cases").joinpath("test_lines.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())
        # set the tolerance to merge the given lines
        geo.epsilon = 1e-6
        # geo.merge_points()
        geo.merge_lines()
        # there is only one described surface exists in the given geometry
        surfaces = geo.find_surfaces()

        # create a gmsh mesh from the given geometry
        gmsh = GMSHModel(geo)
        gmsh.gmsh_writer("test1")

        msh_data = read("test1.msh")

        # tests that the code generates a valid mesh
        self.assertGreaterEqual(len(msh_data.cells[0].data), 3)

        # check the surface, the surface should contain 8 edges
        self.assertEqual(len(surfaces[0]), 8)
        # remove the geo and msh files
        remove("test1.geo_unrolled")
        remove("test1.msh")

    def test_bezier_line_surface(self):
        # import the surface
        eml = files("tests.pygmsh_tests.test_cases").joinpath("test_bezier.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())
        # set the tolerance to merge the given lines
        geo.epsilon = 1e-6
        # geo.merge_points()

        # there is only one described surface exists in the given geometry
        surfaces = geo.find_surfaces()
        geo.plot_connection_graph(debug=True)
        # create a gmsh mesh from the given geometry
        gmsh = GMSHModel(geo)
        gmsh.gmsh_writer("test2")

        msh_data = read("test2.msh")

        # tests that the code generates a valid mesh
        self.assertGreaterEqual(len(msh_data.cells[0].data), 3)

        # check the surface, the surface should contain only 7 edges
        self.assertEqual(len(surfaces[0]), 7)

        # remove the geo and msh files
        remove("test2.geo_unrolled")
        remove("test2.msh")

    def test_circle_defined_surface(self):
        # define the geometry by hand, a simple arc
        geo = Geometry()

        a = Node(x=0.0, y=0.0, id_=1)
        b = Node(x=10.0, y=0.0, id_=2)
        c = Node(x=0.0, y=10.0, id_=3)

        geo.add_line(Line(a, b))
        geo.add_line(Line(a, c))
        geo.add_arc(CircleArc(c, a, b))

        # there is only one described surface exists in the given geometry
        surfaces = geo.find_surfaces()

        gmsh = GMSHModel(geo)
        gmsh.gmsh_writer("test3")

        # check the surface, the surface should contain only 3 edges
        self.assertEqual(len(surfaces[0]), 3)
        self.assertEqual(round(surfaces[0][0].start_pt.x, 1), 0.0)

        # remove the geo and msh files
        remove("test3.geo_unrolled")
        remove("test3.msh")
