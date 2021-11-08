from copy import copy
from math import pi
from unittest import TestCase

from digital_twin_distiller.objects import CircleArc, CubicBezier, Line, Node, ParametricBezier, Rectangle


class TestNodeOperations(TestCase):
    def test_node_init(self):
        # initialize without an id number
        self.assertEqual((1.0, 1.0), Node(1.0, 1.0).as_tuple())

    def test_operators(self):
        n0 = Node(5, -1)
        n1 = Node(-5, 1)
        n0 = n0 + n1
        self.assertEqual(n0, Node(0, 0))
        self.assertEqual(n0 + 6, Node(6, 6))

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
        self.assertAlmostEqual(n0.angle_to(n1), 355 / 113 / 2, 3)


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
        self.assertAlmostEqual(l1.distance_to_point(0.5, 1), 1.0, 3)
        self.assertAlmostEqual(l1.distance_to_point(1, 1), 1.0, 3)
        self.assertAlmostEqual(l1.distance_to_point(2, 0), 1.0, 3)
        self.assertAlmostEqual(l1.distance_to_point(-1, 0), 1.0, 3)


class TestArc(TestCase):
    def test_creation(self):
        c = CircleArc(start_pt=Node(-1, 0), center_pt=Node(0, 0), end_pt=Node(1, 0))

        self.assertEqual(c.start_pt, Node(-1, 0))
        self.assertEqual(c.center_pt, Node(0, 0))
        self.assertEqual(c.end_pt, Node(1, 0))
        self.assertEqual(c.apex_pt, Node(0, -1))
        self.assertAlmostEqual(c.radius, 1.0, 5)

    def test_copy(self):
        c = CircleArc(start_pt=Node(-1, 0), center_pt=Node(0, 0), end_pt=Node(1, 0))
        c1 = copy(c)

        for attr_i in c.__dict__:
            if attr_i != "id":
                self.assertEqual(getattr(c, attr_i), getattr(c1, attr_i))

        self.assertNotEqual(c.id, c1.id)

    def test_distance(self):
        c = CircleArc(start_pt=Node(-1, 0), center_pt=Node(0, 0), end_pt=Node(1, 0))

        self.assertAlmostEqual(c.distance_to_point(0, 0), 1.0, 5)
        self.assertAlmostEqual(c.distance_to_point(5, 0), 4.0, 5)
        self.assertAlmostEqual(c.distance_to_point(-5, 0), 4.0, 5)
        self.assertAlmostEqual(c.distance_to_point(0, -5), 4.0, 5)


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
        bz = ParametricBezier((1, 0), (2, 2), (5, 5), (10, 0))
        self.assertEqual(bz.p0, (1, 0))
        self.assertEqual(bz.p1, (2, 2))
        self.assertEqual(bz.p2, (5, 5))
        self.assertEqual(bz.p3, (10, 0))

    def test_set(self):
        bz = ParametricBezier((1, 0), (2, 2), (5, 5), (10, 0))

        bz.set(start=(0, 0))
        self.assertEqual(bz.p0, (0, 0))
        self.assertEqual(bz.p1, (2, 2))
        self.assertEqual(bz.p2, (5, 5))
        self.assertEqual(bz.p3, (10, 0))

    def test_approximate(self):
        bz = ParametricBezier((0, 0), (2.5, 0), (7.5, 0), (10, 0))

        reflines = [
            Line(Node(0, 0), Node(5.0, 0)),
            Line(Node(5.0, 0), Node(10.0, 0)),
        ]
        for refl, l in zip(reflines, bz.approximate(2)):
            self.assertEqual(refl, l)

    def test_call(self):
        bz = ParametricBezier((0, 0), (2.5, 0), (7.5, 0), (10, 0))
        x, y = bz(0)
        self.assertAlmostEqual(x, 0.0, 5)
        self.assertAlmostEqual(y, 0.0, 5)

        x, y = bz(0.5)
        self.assertAlmostEqual(x, 5.0, 5)
        self.assertAlmostEqual(y, 0.0, 5)

        x, y = bz(1)
        self.assertAlmostEqual(x, 10.0, 5)
        self.assertAlmostEqual(y, 0.0, 5)

        with self.assertRaises(AssertionError):
            x, y = bz(5)


class TesRectangle(TestCase):
    def test_creation(self):
        r = Rectangle(x0=0, y0=0, width=1, height=2)
        self.assertEqual(r.a, Node(0, 0))
        self.assertEqual(r.b, Node(1, 0))
        self.assertEqual(r.c, Node(1, 2))
        self.assertEqual(r.d, Node(0, 2))

        r = Rectangle(x0=0, y0=0, x1=1, y1=2)
        self.assertEqual(r.a, Node(0, 0))
        self.assertEqual(r.b, Node(1, 0))
        self.assertEqual(r.c, Node(1, 2))
        self.assertEqual(r.d, Node(0, 2))

        with self.assertRaises(ValueError):
            r = Rectangle(x0=0, y0=0, width=5)

    def test_rotation(self):
        r = Rectangle(x0=0, y0=0, width=1, height=2)
        r._print_sidelengths()
        r.rotate(360)
        r.rotate(180, fx_point="a")
        r.rotate(180, fx_point="b")
        r.rotate(180, fx_point="c")
        r.rotate(180, fx_point="d")
        print(r)

        self.assertEqual(r.a, Node(0, 0))
        self.assertEqual(r.b, Node(1, 0))
        self.assertEqual(r.c, Node(1, 2))
        self.assertEqual(r.d, Node(0, 2))

        with self.assertRaises(ValueError):
            r.rotate(3, fx_point="falspoint")

    def test_set_width(self):
        r = Rectangle(x0=0, y0=0, width=1, height=2)

        r.set_width(5, fx_point="a")
        r.set_width(10, fx_point="c")
        r.set_width(1)

        self.assertEqual(r.a, Node(-0.5, 0))
        self.assertEqual(r.b, Node(0.5, 0))
        self.assertEqual(r.c, Node(0.5, 2))
        self.assertEqual(r.d, Node(-0.5, 2))

        with self.assertRaises(ValueError):
            r.set_width(1, fx_point="false_point")

    def test_set_height(self):
        r = Rectangle(x0=0, y0=0, width=1, height=2)

        r.set_height(5, fx_point="a")
        r.set_height(10, fx_point="c")
        r.set_height(1)

        self.assertEqual(r.a, Node(0, -0.5))
        self.assertEqual(r.b, Node(1, -0.5))
        self.assertEqual(r.c, Node(1, 0.5))
        self.assertEqual(r.d, Node(0, 0.5))

        with self.assertRaises(ValueError):
            r.set_height(1, fx_point="false_point")

    def test_put(self):
        r = Rectangle(x0=0, y0=0, width=1, height=2)

        r.put(10, 40)
        r.put(40, -20, fx_point="b")
        r.put(7, 0.76, fx_point="c")
        r.put(-65, 5555, fx_point="d")
        r.put(0, 0, fx_point="a")

        self.assertEqual(r.a, Node(0, 0))
        self.assertEqual(r.b, Node(1, 0))
        self.assertEqual(r.c, Node(1, 2))
        self.assertEqual(r.d, Node(0, 2))

        with self.assertRaises(ValueError):
            r.put(0, 0, fx_point="invalid")

    def test_mirror(self):
        r = Rectangle(x0=1, y0=1, width=1, height=2)
        r.mirror(p1=(0, 0), p2=(0, 1))

        self.assertEqual(r.a, Node(-1, 1))
        self.assertEqual(r.b, Node(-2, 1))
        self.assertEqual(r.c, Node(-2, 3))
        self.assertEqual(r.d, Node(-1, 3))

    def test_copy(self):
        r = Rectangle(x0=0, y0=0, width=1, height=2)
        r1 = copy(r)

        self.assertEqual(r.a, r1.a)
        self.assertEqual(r.b, r1.b)
        self.assertEqual(r.c, r1.c)
        self.assertEqual(r.d, r1.d)
        self.assertAlmostEqual(r.width, r1.width, 5)
        self.assertAlmostEqual(r.height, r1.height, 5)


if __name__ == "__main__":
    t = TestParametricBezier()
    t.test_call()
