import os
from unittest import TestCase

from adze_modeler.agros2d_wrapper import Agros2DWrapper


class FemmTester(TestCase):
    def get_fresh_wrapper(self):
        return Agros2DWrapper()

    def test_set_problem_type(self):
        # TODO: test for other types of problems
        w = self.get_fresh_wrapper()

        w.set_problem_type("electrostatic")
        self.assertEqual(w.problem_type, "electrostatic")

        w.set_problem_type("magnetic")
        self.assertEqual(w.problem_type, "magnetic")

        # test for nonexisting field
        self.assertRaises(NotImplementedError, w.set_problem_type, "eper")
