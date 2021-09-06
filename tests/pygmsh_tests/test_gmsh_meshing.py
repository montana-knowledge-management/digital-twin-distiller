from unittest import TestCase
from adze_modeler.geometry import Geometry
from adze_modeler.gmsh import GMSHModel

from importlib_resources import files
from os import remove
from meshio._helpers import read


# # plotting out the mesh
# import pyvista as pv

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

        # there is only one described surface exists in the given geometry
        surfaces = geo.find_surfaces()
        geo.export_svg('test_cases/test_lines.svg')
        # create a gmsh mesh from the given geometry
        gmsh = GMSHModel(geo)
        gmsh.gmsh_writer('test1')

        msh_data = read('test1.vtk')

        # tests that the code generates a valid mesh
        self.assertGreaterEqual(len(msh_data.cells), 3)

        # check the surface, the surface should contain 8 edges
        self.assertEqual(len(surfaces[0]),8)
        self.assertEqual(round(surfaces[0][0].start_pt.x,1), 103.4)
        # remove the geo and msh files
        remove('test1.geo_unrolled')
        remove('test1.vtk')
