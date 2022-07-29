from unittest import TestCase

from digital_twin_distiller.geometry import Geometry

test_dict = {
    "style": "fill:none;stroke:#ff0000;stroke-width:0.265;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;stroke-miterlimit:4;stroke-dasharray:none",
    "d": "M 92.817983,139.46116 C 125.27176,102.78228 137.94025,122.3735 147.06331,152.04803",
    "id": "path10",
    "inkscape:connector-curvature": "0",
    "sodipodi:nodetypes": "cc",
}


class TestColor(TestCase):
    def test_example_dict(self):
        color = Geometry.get_color_value_from_svg(test_dict)
        self.assertEqual(color, "#ff0000")

    def test_no_color(self):
        color = Geometry.get_color_value_from_svg({})
        self.assertEqual(color, "#000000")
