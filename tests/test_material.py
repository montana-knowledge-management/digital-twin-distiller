import unittest
from adze_modeler.material import Material
from copy import copy


class TestMaterial(unittest.TestCase):
    def test_material(self):
        m = Material("air")
        self.assertEqual(m.name, "air")

        m1 = copy(m)
        self.assertEqual(m1.name, m.name)
