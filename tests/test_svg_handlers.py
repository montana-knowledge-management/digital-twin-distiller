from unittest import TestCase

from adze_modeler.geometry import Geometry
from importlib_resources import files


class TestSvgImport(TestCase):
    def test_owl_import_to_geometry(self):
        eml = files("examples.owl").joinpath("owl-svgrepo-com.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())

        # checks the first coordinate of the first node
        self.assertEqual(445.642, geo.nodes[0].x)
        # self.assertEqual(635, geo.nodes[-1].id)
        # the number of lines and cubicbeziers should be larger than 0
        self.assertTrue(len(geo.lines) > 0)
        self.assertTrue(len(geo.cubic_beziers) > 0)
        self.assertTrue(len(geo.circle_arcs) == 0)

    def test_approximate_owl(self):
        eml = files("examples.owl").joinpath("owl-svgrepo-com.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())

        print(geo.cubic_beziers)

