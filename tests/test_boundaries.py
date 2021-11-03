import unittest

from adze_modeler.boundaries import (
    DirichletBoundaryCondition,
    NeumannBoundaryCondition,
)


class TestDirichletBoundary(unittest.TestCase):
    def get_bc(self, name, field_type):
        return DirichletBoundaryCondition(name=name, field_type=field_type)

    def test_proper_field_type(self):

        bc = self.get_bc("alma", "magnetic")
        self.assertEqual("magnetic", bc.field)

        self.assertRaises(
            ValueError,
            DirichletBoundaryCondition,
            name="alma",
            field_type="DummyFieldType",
        )

    def test_proper_representation(self):
        bc = self.get_bc("alma", "magnetic")
        bc.set_value("magnetic_potential", 2)
        self.assertEqual(
            bc.__str__(),
            "name: alma, type: magnetic-dirichlet, value(s): magnetic_potential: 2",
        )

    def test_init_with_values(self):
        bc = DirichletBoundaryCondition(
            "alma", "electrostatic", fixed_voltage=10
        )
        self.assertEqual(bc.valuedict["fixed_voltage"], 10)

        self.assertRaises(
            ValueError,
            DirichletBoundaryCondition,
            name="alma",
            field_type="heat",
            magnetic_potential=10,
        )

    def test_magnetic_dirichlet_bc(self):
        bc = self.get_bc("a0", "magnetic")
        self.assertTrue("magnetic_potential" in bc.accepted_keys["magnetic"])

        bc.set_value("magnetic_potential", 42)
        self.assertEqual(42, bc.valuedict["magnetic_potential"])

        self.assertRaises(ValueError, bc.set_value, key="falsekey", value=-1)

    def test_electrostatic_dirichlet_bc(self):
        bc = self.get_bc("Vg", "electrostatic")
        self.assertTrue("fixed_voltage" in bc.accepted_keys["electrostatic"])

        bc.set_value("fixed_voltage", 42)
        self.assertEqual(42, bc.valuedict["fixed_voltage"])

    def test_current_dirichlet_bc(self):
        bc = self.get_bc("Vg", "current")
        self.assertTrue("fixed_voltage" in bc.accepted_keys["electrostatic"])

        bc.set_value("fixed_voltage", 42)
        self.assertEqual(42, bc.valuedict["fixed_voltage"])

    def test_heat_dirichlet_bc(self):
        bc = self.get_bc("Tf", "heat")
        self.assertTrue("temperature" in bc.accepted_keys["heat"])

        bc.set_value("temperature", 42)
        self.assertEqual(42, bc.valuedict["temperature"])


class TestNeumannBoundary(unittest.TestCase):
    def get_bc(self, name, field_type):
        return NeumannBoundaryCondition(name=name, field_type=field_type)

    def test_proper_field_type(self):

        bc = self.get_bc("alma", "magnetic")
        self.assertEqual("magnetic", bc.field)

        self.assertRaises(
            ValueError,
            NeumannBoundaryCondition,
            name="alma",
            field_type="DummyFieldType",
        )

    def test_proper_representation(self):
        bc = self.get_bc("alma", "magnetic")
        bc.set_value("surface_current", 2)
        self.assertEqual(
            bc.__str__(),
            "name: alma, type: magnetic-neumann, value(s): surface_current: 2",
        )

    def test_init_with_values(self):
        bc = NeumannBoundaryCondition(
            "alma", "electrostatic", surface_charge_density=10
        )
        self.assertEqual(bc.valuedict["surface_charge_density"], 10)

        self.assertRaises(
            ValueError,
            NeumannBoundaryCondition,
            name="alma",
            field_type="heat",
            magnetic_potential=10,
        )

    def test_magnetic_neumann_bc(self):
        bc = self.get_bc("a0", "magnetic")
        self.assertTrue("surface_current" in bc.accepted_keys["magnetic"])

        bc.set_value("surface_current", 42)
        self.assertEqual(42, bc.valuedict["surface_current"])

        self.assertRaises(ValueError, bc.set_value, key="falsekey", value=-1)

    def test_electrostatic_neumann_bc(self):
        bc = self.get_bc("Vg", "electrostatic")
        self.assertTrue(
            "surface_charge_density" in bc.accepted_keys["electrostatic"]
        )

        bc.set_value("surface_charge_density", 42)
        self.assertEqual(42, bc.valuedict["surface_charge_density"])

    def test_current_neumann_bc(self):
        bc = self.get_bc("Vg", "current")
        self.assertTrue("current_density" in bc.accepted_keys["current"])

        bc.set_value("current_density", 42)
        self.assertEqual(42, bc.valuedict["current_density"])

    def test_heat_neumann_bc(self):
        bc = self.get_bc("Tf", "heat")
        self.assertTrue("heat_flux" in bc.accepted_keys["heat"])

        bc.set_value("heat_flux", 42)
        self.assertEqual(42, bc.valuedict["heat_flux"])
