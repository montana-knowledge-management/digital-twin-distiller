from pathlib import Path
from unittest import TestCase

from digital_twin_distiller.geometry import Geometry

resources = Path(__file__).parent


class TestSvgImport(TestCase):
    def test_owl_import_to_geometry(self):
        eml = Path(resources / "svgtests" / "owl-shape.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())

        # checks the first coordinate of the first node
        # self.assertTrue(abs(442.466 - geo.nodes[0].x) < 1e-3)
        # self.assertEqual(635, geo.nodes[-1].id)
        # the number of lines and cubicbeziers should be larger than 0
        self.assertTrue(len(geo.lines) > 0)
        self.assertTrue(len(geo.cubic_beziers) > 0)
        self.assertTrue(len(geo.circle_arcs) == 0)

    def test_approximate_owl(self):
        eml = Path(resources / "svgtests" / "owl-shape.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())

        print(geo.cubic_beziers)
