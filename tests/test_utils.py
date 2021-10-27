import unittest

import adze_modeler.utils as u
from adze_modeler.objects import Node


class TestUtils(unittest.TestCase):
    def test_get_id(self):
        self.assertTrue(u.getID() > 0)

    def test_mirror(self):
        p0 = Node(0, 0)
        p1 = Node(0, 1)
        p2 = Node(-1, 1)

        p3 = u.mirror_point(p0, p1, p2)
        self.assertEqual(p3, Node(1, 1))
