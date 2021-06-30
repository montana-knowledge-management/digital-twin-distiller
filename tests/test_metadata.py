from adze_modeler.metadata import Agros2DMetadata
from adze_modeler.metadata import FemmMetadata

import unittest

class TestAgros2DMetadata(unittest.TestCase):
    def test_metadata(self):

        m = Agros2DMetadata()
        self.assertEqual(m.file_suffix, ".py")
        self.assertRaises(SystemExit, m.validate_metadata)

        m = FemmMetadata()
        m.file_script_name="test.fec"
        self.assertEqual(m.file_suffix, ".lua")
        self.assertEqual(m.validate_metadata(), None)
        m.unit = 1e-3
        self.assertRaises(ValueError, m.validate_metadata)