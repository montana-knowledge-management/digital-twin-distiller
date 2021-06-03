from math import pi
from unittest import TestCase

from adze_modeler.objects import CubicBezier
from adze_modeler.objects import Line
from adze_modeler.objects import Node


class TestNodeOperations(TestCase):
    def test_node_init(self):
        # initialize without an id number
        self.assertEqual((1.0, 1.0), Node(1.0, 1.0).as_tuple())

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
        self.assertEqual("(1.0, 0.0, id=None,label=None)", str(a))

    def test_distance(self):
        a = Node(1.0, 0.0)
        b = Node(0.5, 0.0)
        c = a + b
        d = a - b
        e = a * 2.0
        self.assertEqual((1.5, 0.0), c.as_tuple())
        self.assertEqual((0.5, 0.0), d.as_tuple())
        self.assertEqual((2.0, 0.0), e.as_tuple())


class TestLine(TestCase):
    def test_init_line(self):
        a = Node(1.0, 0.0)
        b = Node(0.5, 0.0)

        l = Line(a, b, id=1, label="test")

        self.assertEqual("test", l.label)
        self.assertEqual(a, l.start_pt)
        self.assertEqual(b, l.end_pt)

        # repr string
        self.assertIn(
            "Line(Node(1.0, 0.0, id=None,label=None), Node(0.5, 0.0, id=None,label=None), id=1,label='test')", str(l)
        )


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
            "CubicBezier(Node(1.0, 0.0, id=None,label=None), Node(0.6, 0.1, id=None,label=None), Node("
            + "0.7, 0.2, id=None,label=None), Node(0.5, 0.0, id=None,label=None), id=1,label='test')",
            str(cb),
        )
        # print(cb)
