import adze_modeler.utils as u
import unittest


class TestUtils(unittest.TestCase):
    def test_get_id(self):
        self.assertTrue(u.getID() > 0)
