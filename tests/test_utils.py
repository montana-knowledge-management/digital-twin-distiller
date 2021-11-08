import unittest
import pytest
import os

import digital_twin_distiller.utils as u
from digital_twin_distiller.objects import Node


class TestUtils(unittest.TestCase):
    def test_get_id(self):
        self.assertTrue(u.getID() > 0)

    def test_mirror(self):
        p0 = Node(0, 0)
        p1 = Node(0, 1)
        p2 = Node(-1, 1)

        p3 = u.mirror_point(p0, p1, p2)
        self.assertEqual(p3, Node(1, 1))

    def test_mm2px(self):
        input1 = 0
        input2 = 1
        input3 = 2.75
        # negative number?

        coverted1 = u.mm2px(input1)
        coverted2 = u.mm2px(input2)
        coverted3 = u.mm2px(input3)

        self.assertEqual(coverted1, 0)
        self.assertEqual(coverted2, int(3.7795275591))
        self.assertEqual(coverted3, int(10.393700787))

    def test_mm2inch(self):
        input1 = 0
        input2 = 1
        input3 = 25.4

        coverted1 = u.mm2inch(input1)
        coverted2 = u.mm2inch(input2)
        coverted3 = u.mm2inch(input3)

        self.assertEqual(coverted1, 0)
        self.assertEqual(coverted2, 0.03937007874)
        self.assertEqual(coverted3, 0.999999999996)

    def test_get_width_height(self):
        with pytest.raises(ValueError):
            u.get_width_height(type_="notexists", aspect=(16, 10), unit="px")

        result = u.get_width_height(
            type_="minimal",
            aspect=(16, 9),
            unit="mm"
        )
        scale = 16 / 9
        self.assertEqual(result, (30, 30 / scale))  # 113, 63

    def test_inch2mm(self):
        result = u.inch2mm(1)

        self.assertEqual(25.4, result)

    def test_rms(self):
        result = u.rms((1, 1, 1))
        self.assertEqual(result, 1)

    def test_pairwise(self):
        case1 = ('A')

        self.assertTrue(u.pairwise(case1))
        self.assertEqual(list(u.pairwise(case1)), [])

        case2 = '1234'
        self.assertEqual(list(u.pairwise(case2)), [('1', '2'), ('2', '3'), ('3', '4')])

    def test_csv_read(self):
        my_path = os.path.abspath(os.path.dirname(__file__))
        path = os.path.join(my_path, "test_utils/test1.csv")

        test = u.csv_read(path, True)
        self.assertEqual(test, {1.0: (245.0, 54.0, 54.0), 2.0: (12.0, 54.0, 123.0)})

    def test_polyfit_line_to_line(self):
        x = [1, 2, 3]
        y = [1, 2, 3]

        x_fine, y_fine = u.get_polyfit(x, y)

        self.assertGreater(len(x_fine), len(x))
        self.assertEqual(round(x_fine[1], 3), 1.333)
        self.assertEqual(round(y_fine[1], 3), 1.764)
