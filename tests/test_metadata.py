import unittest
from copy import copy

from adze_modeler.metadata import Agros2DMetadata
from adze_modeler.metadata import FemmMetadata


class TestMetadata(unittest.TestCase):
    def test_metadata(self):

        m = Agros2DMetadata()
        self.assertEqual(m.file_suffix, ".py")
        self.assertRaises(SystemExit, m.validate_metadata)

        m = FemmMetadata()
        m.file_script_name = "test.fec"
        self.assertEqual(m.file_suffix, ".lua")
        self.assertEqual(m.validate_metadata(), None)
        m.unit = 1e-3
        self.assertRaises(ValueError, m.validate_metadata)

    def test_copy(self):

        agros_metadata = Agros2DMetadata()
        agros_metadata.file_script_name = "test1"
        agros_metadata.file_metrics_name = "test2"
        agros_metadata.problem_type = "magnetic"
        agros_metadata.coordinate_type = "axisymmetric"
        agros_metadata.analysis_type = "steadystate"
        agros_metadata.unit = 1e-3
        agros_metadata.nb_refinements = 0
        agros_metadata.adaptivity = "hp-adaptivity"
        agros_metadata.adaptivity_tol = 1

        metadata1 = copy(agros_metadata)

        for attr_i in metadata1.__dict__:
            self.assertEqual(getattr(agros_metadata, attr_i), getattr(metadata1, attr_i))

        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "axisymmetric"
        femm_metadata.file_script_name = "alma"
        femm_metadata.file_metrics_name = "banan"
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False

        metadata2 = copy(femm_metadata)

        for attr_i in metadata2.__dict__:
            self.assertEqual(getattr(femm_metadata, attr_i), getattr(metadata2, attr_i))
