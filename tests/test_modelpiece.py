import unittest

from adze_modeler.modelpiece import ModelPiece
from adze_modeler.objects import CircleArc, Line, Node


class TestModelpiece(unittest.TestCase):
    def get_modelpiece(self):
        m = ModelPiece('test')
        n0 = Node(0, 1)
        l0 = Line(Node(0, 0), Node(0, -1))
        c0 = CircleArc(start_pt=Node(-1, 0),
                       center_pt=Node(0, 0),
                       end_pt=Node(1, 0))

        self.assertEqual(c0.apex_pt, Node(0, -1))

        m.geom.add_node(n0)
        m.geom.add_line(l0)
        m.geom.add_arc(c0)
        m.update_bbox()
        return m

    def test_bbox(self):
        m = self.get_modelpiece()

        self.assertAlmostEqual(m.bbox[0], -1.0, delta=1e-9)
        self.assertAlmostEqual(m.bbox[1], -1.0, delta=1e-9)
        self.assertAlmostEqual(m.bbox[2], 1.0, delta=1e-9)
        self.assertAlmostEqual(m.bbox[3], 1.0, delta=1e-9)

    def test_translate(self):
        m = self.get_modelpiece()

        m.translate(5, -1)

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(5, 0))
        self.assertEqual(m.geom.lines[0].start_pt, Node(5, -1))
        self.assertEqual(m.geom.lines[0].end_pt, Node(5, -2))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(4, -1))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(5, -1))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(5, -2))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(6, -1))

    def test_put_lower_left(self):
        m = self.get_modelpiece()


        m.put(10, 10, bbox_ref="lower-left")

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(11, 12))
        self.assertEqual(m.geom.lines[0].start_pt, Node(11, 11))
        self.assertEqual(m.geom.lines[0].end_pt, Node(11, 10))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(10, 11))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(11, 11))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(11, 10))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(12, 11))

    def test_put_lower_right(self):
        m = self.get_modelpiece()


        m.put(10, 10, bbox_ref="lower-right")

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(9, 12))
        self.assertEqual(m.geom.lines[0].start_pt, Node(9, 11))
        self.assertEqual(m.geom.lines[0].end_pt, Node(9, 10))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(8, 11))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(9, 11))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(9, 10))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(10, 11))

    def test_put_upper_left(self):
        m = self.get_modelpiece()


        m.put(10, 10, bbox_ref="upper-left")

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(11, 10))
        self.assertEqual(m.geom.lines[0].start_pt, Node(11, 9))
        self.assertEqual(m.geom.lines[0].end_pt, Node(11, 8))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(10, 9))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(11, 9))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(11, 8))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(12, 9))

    def test_put_upper_right(self):
        m = self.get_modelpiece()

        m.put(10, 10, bbox_ref="upper-right")

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(9, 10))
        self.assertEqual(m.geom.lines[0].start_pt, Node(9, 9))
        self.assertEqual(m.geom.lines[0].end_pt, Node(9, 8))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(8, 9))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(9, 9))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(9, 8))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(10, 9))

    def test_put_upper(self):
        m = self.get_modelpiece()

        m.put(10, 10, bbox_ref="upper")

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(10, 10))
        self.assertEqual(m.geom.lines[0].start_pt, Node(10, 9))
        self.assertEqual(m.geom.lines[0].end_pt, Node(10, 8))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(9, 9))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(10, 9))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(10, 8))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(11, 9))

    def test_put_lower(self):
        m = self.get_modelpiece()

        m.put(10, 10, bbox_ref="lower")

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(10, 12))
        self.assertEqual(m.geom.lines[0].start_pt, Node(10, 11))
        self.assertEqual(m.geom.lines[0].end_pt, Node(10, 10))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(9, 11))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(10, 11))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(10, 10))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(11, 11))

    def test_put_right(self):
        m = self.get_modelpiece()

        m.put(10, 10, bbox_ref="right")

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(9, 11))
        self.assertEqual(m.geom.lines[0].start_pt, Node(9, 10))
        self.assertEqual(m.geom.lines[0].end_pt, Node(9, 9))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(8, 10))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(9, 10))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(9, 9))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(10, 10))

    def test_put_left(self):
        m = self.get_modelpiece()

        m.put(10, 10, bbox_ref="left")

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(11, 11))
        self.assertEqual(m.geom.lines[0].start_pt, Node(11, 10))
        self.assertEqual(m.geom.lines[0].end_pt, Node(11, 9))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(10, 10))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(11, 10))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(11, 9))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(12, 10))

    def test_rotate(self):
        m = self.get_modelpiece()

        m.rotate(ref_point=(0, 0), alpha=90)

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(-1, 0))
        self.assertEqual(m.geom.lines[0].start_pt, Node(0, 0))
        self.assertEqual(m.geom.lines[0].end_pt, Node(1, 0))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(0, -1))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(0, 0))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(1, 0))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(0, 1))

    def test_scale(self):
        m = self.get_modelpiece()

        m.scale(10, 10)

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(0, 10))
        self.assertEqual(m.geom.lines[0].start_pt, Node(0, 0))
        self.assertEqual(m.geom.lines[0].end_pt, Node(0, -10))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(-10, 0))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(0, 0))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(0, -10))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(10, 0))

    def test_mirror(self):
        m = self.get_modelpiece()

        m.put(0, 10, bbox_ref="lower")
        m.mirror((0, 0), (1, 0))

        self.assertEqual(len(m.geom.nodes), 5)
        self.assertEqual(m.geom.nodes[0], Node(0, -12))
        self.assertEqual(m.geom.lines[0].start_pt, Node(0, -11))
        self.assertEqual(m.geom.lines[0].end_pt, Node(0, -10))
        self.assertEqual(m.geom.circle_arcs[0].start_pt, Node(1, -11))
        self.assertEqual(m.geom.circle_arcs[0].center_pt, Node(0, -11))
        self.assertEqual(m.geom.circle_arcs[0].apex_pt, Node(0, -10))
        self.assertEqual(m.geom.circle_arcs[0].end_pt, Node(-1, -11))
