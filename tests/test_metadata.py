import unittest
from copy import copy

from digital_twin_distiller.metadata import Agros2DMetadata, FemmMetadata


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

    def test_agros_should_pass_file(self):
        self.should_pass_file(Agros2DMetadata())

    def test_femm_should_pass_file(self):
        self.should_pass_file(FemmMetadata())

    def should_pass_file(self, m):

        path_prefix = "t/t2.t3/t4/t5@t6/t7"

        self.file_script_name_assertion(m, path_prefix, "name.txt")
        self.file_script_name_assertion(m, path_prefix, "name")
        self.file_script_name_assertion(m, path_prefix, "name.test.txt")

        path_prefix = "t/t2/t3/t4"

        self.file_script_name_assertion(m, path_prefix, "name.txt")
        self.file_script_name_assertion(m, path_prefix, "name")
        self.file_script_name_assertion(m, path_prefix, "name.test.txt")

    def file_script_name_assertion(self, m, path_prefix, file):

        m.file_script_name = path_prefix + "/" + file

        m.validate_file_name()
        self.assertNotEqual(file, m.file_script_name)

        file_name = path_prefix + "/" + (file if file.rfind(".") == -1 else file[: file.rfind(".")])

        self.assertEqual(file_name + m.file_suffix, m.file_script_name)
