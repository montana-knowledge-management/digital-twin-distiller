from unittest import TestCase

from adze_modeler.geometry import Geometry
from importlib_resources import files


class TestDXFImport(TestCase):
    def test_dxf_import_to_geometry(self):
        eml = files("examples.induction_motor").joinpath("2horse.dxf")
        # print(eml)

        geo = Geometry()
        geo.import_dxf(eml.as_posix())

        self.assertEqual(236, len(geo.nodes))
        self.assertEqual(89, len(geo.lines))
