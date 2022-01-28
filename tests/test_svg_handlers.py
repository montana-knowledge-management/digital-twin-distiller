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
        # the number of lines and cubicbeziers should be larger than 0
        self.assertTrue(len(geo.lines) > 0)
        self.assertTrue(len(geo.cubic_beziers) > 0)
        self.assertTrue(len(geo.circle_arcs) == 0)

    def test_approximate_owl(self):
        eml = Path(resources / "svgtests" / "owl-shape.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())

        # print(geo.cubic_beziers)
        self.assertEqual(len(geo.cubic_beziers), 46)

    def test_bez_approx(self):
        eml = Path(resources / "svgtests" / "test_bez_approx.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())

        # the line should be red
        self.assertEqual(geo.lines[0].color, '#ff00ff')
        self.assertEqual(geo.cubic_beziers[0].color, '#ff0000')

    def test_colored_rotor_import(self):
        eml = Path(resources / "svgtests" / "antunes_rotor.svg")
        geo = Geometry()
        geo.import_svg(eml.as_posix())

        colors = []
        for c in geo.circle_arcs:
            colors.append(c.color)

        # yellow
        self.assertIn('#ffcc00', colors)
        # blue
        self.assertIn('#0000ff', colors)
        # black
        self.assertIn('#000000', colors)
        # aqua
        self.assertIn('#00ff00', colors)

        line_colors = []
        for lc in geo.lines:
            line_colors.append(line_colors.append(lc.color))
        # red
        self.assertIn('#ff0000', line_colors)
