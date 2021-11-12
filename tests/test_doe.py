import unittest

from digital_twin_distiller.doe import doe_bbdesign, doe_ccf, doe_fullfact, doe_pbdesign


class TestDOE(unittest.TestCase):
    def test_boxbehnkendesign(self):
        self.assertEqual(doe_bbdesign(3, center=1)[0], [-1, -1, 0])
        self.assertEqual(len(doe_bbdesign(3, center=1)), 13)

    def test_fullfactorial(self):
        self.assertEqual(len(doe_fullfact([3] * 4)), 81)

    def test_pbdesign(self):
        self.assertEqual(len(doe_pbdesign(3)), 4)

    def test_ccf(self):
        self.assertEqual(len(doe_ccf(4)), 25)
