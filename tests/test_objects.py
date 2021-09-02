from copy import copy

from adze_modeler.objects import CubicBezier, ParametricBezier
from adze_modeler.objects import Line
from adze_modeler.objects import Node
from math import pi
from unittest import TestCase


class TestNodeOperations(TestCase):
    def test_node_init(self):
        # initialize without an id number
        self.assertEqual((1.0, 1.0), Node(1.0, 1.0).as_tuple())

    def test_operators(self):
        n0 = Node(5, -1)
        n1 = Node(-5, 1)
        n0 = n0 + n1
        self.assertEqual(n0, Node(0, 0))
        self.assertEqual(n0+6, Node(6, 6))

    def test_iter(self):
        n0 = Node(9, 6)
        l = list(n0)
        self.assertAlmostEqual(l[0], n0.x, 5)
        self.assertAlmostEqual(l[1], n0.y, 5)

    def test_rotate_a_node(self):
        a = Node(1.0, 0.0)
        c = a.rotate(pi / 2)
        self.assertEqual((0.0, 1.0), c.as_tuple())

    def test_rotate_node_around_a_point(self):
        a = Node(1.0, 0.0)
        b = Node(0.5, 0.0)

        d = a.rotate_about(b, pi / 2)
        self.assertEqual((0.5, 0.5), d.as_tuple())

    def test_strings(self):
        a = Node(1.0, 0.0)
        # self.assertEqual(f"(1.0, 0.0, label=None)", str(a))

    def test_distance(self):
        a = Node(1.0, 0.0)
        b = Node(0.5, 0.0)
        c = a + b
        d = a - b
        e = a * 2.0
        self.assertEqual((1.5, 0.0), c.as_tuple())
        self.assertEqual((0.5, 0.0), d.as_tuple())
        self.assertEqual((2.0, 0.0), e.as_tuple())

    def test_unitvector(self):
        n0 = Node(0, 0)
        n1 = Node(1, 0)
        u = n0.unit_to(n1)
        self.assertAlmostEqual(u.x, 1.0, 5)
        self.assertAlmostEqual(u.y, 0.0, 5)
        self.assertAlmostEqual(abs(u), 1.0, 5)

    def test_angles(self):
        n0 = Node(1, 0)
        n1 = Node(0, 1)
        self.assertAlmostEqual(n0.angle_to(n1), 355/113/2, 3)


class TestLine(TestCase):
    def test_init_line(self):
        a = Node(1.0, 0.0)
        b = Node(0.5, 0.0)

        l = Line(a, b, id=1, label="test")

        self.assertEqual("test", l.label)
        self.assertEqual(a, l.start_pt)
        self.assertEqual(b, l.end_pt)

        # repr string
        # self.assertIn(f"Line((1.0, 0.0, label=None), (0.5, 0.0, label=None),label='test')", str(l))

    def test_copy_lines(self):
        l1 = Line(Node(0, 0), Node(10, 10))
        l2 = copy(l1)
        self.assertEqual(l1, l2)
        self.assertFalse(l1.start_pt is l2.start_pt)
        self.assertFalse(l1.end_pt is l2.end_pt)

    def test_distance_between_lines(self):
        l1 = Line(Node(0, 0), Node(1, 0))
        self.assertAlmostEqual(l1.distance_to_point(0, 1), 1.0, 3)


class TestCubicBezier(TestCase):
    def test_init_bezier(self):
        a = Node(1.0, 0.0)
        b = Node(0.5, 0.0)

        c1 = Node(0.6, 0.1)
        c2 = Node(0.7, 0.2)

        cb = CubicBezier(a, c1, c2, b, id=1, label="test")

        self.assertEqual("test", cb.label)
        self.assertEqual(a, cb.start_pt)
        self.assertEqual(b, cb.end_pt)

        self.assertIn(
            f"CubicBezier(Node(1.0, 0.0, id={a.id},label=None), Node(0.6, 0.1, id={c1.id},label=None), Node("
            + f"0.7, 0.2, id={c2.id},label=None), Node(0.5, 0.0, id={b.id},label=None), id=1,label='test')",
            str(cb),
        )
        # print(cb)


class TestParametricBezier(TestCase):
    def test_init(self):
        bz = ParametricBezier(start_pt=(1, 0),
                              end_pt=(10, 0),
                              c1=(2, 2),
                              c2=(5, 5))
        self.assertEqual(bz.p0, (1, 0))
        self.assertEqual(bz.p1, (2, 2))
        self.assertEqual(bz.p2, (5, 5))
        self.assertEqual(bz.p3, (10, 0))

    def test_set(self):
        bz = ParametricBezier(start_pt=(1, 0),
                              end_pt=(10, 0),
                              c1=(2, 2),
                              c2=(5, 5))

        bz.set(start_pt=(0, 0))
        self.assertEqual(bz.p0, (0, 0))
        self.assertEqual(bz.p1, (2, 2))
        self.assertEqual(bz.p2, (5, 5))
        self.assertEqual(bz.p3, (10, 0))
