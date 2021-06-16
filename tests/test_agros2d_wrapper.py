import os
from unittest import TestCase

from adze_modeler.agros2d_wrapper import Agros2DWrapper


class Agros2DTester(TestCase):

    def test_add_field_electrostatic(self):
        w = Agros2DWrapper()

        self.assertRaises(ValueError, w.add_field_electrostatic, analysis="eper")
        self.assertRaises(ValueError, w.add_field_electrostatic, solver="eper")
        self.assertRaises(ValueError, w.add_field_electrostatic, matrix_solver="eper")
        self.assertRaises(ValueError, w.add_field_electrostatic, nb_refinements=0.2)
        self.assertRaises(ValueError, w.add_field_electrostatic, nb_polyorder=-4)
        self.assertRaises(ValueError, w.add_field_electrostatic, adaptivity="eper")

        w.add_field_electrostatic(nb_polyorder=3)
        self.assertEqual(w.fields['electrostatic']['nb_polyorder'], 3)

    def test_add_field_magnetic(self):
        w = Agros2DWrapper()

        self.assertRaises(ValueError, w.add_field_magnetic, analysis="eper")
        self.assertRaises(ValueError, w.add_field_magnetic, solver="eper")
        self.assertRaises(ValueError, w.add_field_magnetic, matrix_solver="eper")
        self.assertRaises(ValueError, w.add_field_magnetic, nb_refinements=0.2)
        self.assertRaises(ValueError, w.add_field_magnetic, nb_polyorder=-4)
        self.assertRaises(ValueError, w.add_field_magnetic, adaptivity="eper")

        w.add_field_magnetic(analysis='transient')
        self.assertEqual(w.fields['magnetic']['analysis'], 'transient')

    def test_add_field_heat_transfer(self):
        w = Agros2DWrapper()

        self.assertRaises(ValueError, w.add_field_heat_transfer, analysis="eper")
        self.assertRaises(ValueError, w.add_field_heat_transfer, solver="eper")
        self.assertRaises(ValueError, w.add_field_heat_transfer, matrix_solver="eper")
        self.assertRaises(ValueError, w.add_field_heat_transfer, nb_refinements=0.2)
        self.assertRaises(ValueError, w.add_field_heat_transfer, nb_polyorder=-4)
        self.assertRaises(ValueError, w.add_field_heat_transfer, adaptivity="eper")

        w.add_field_heat_transfer(analysis='transient')
        self.assertEqual(w.fields['heat']['analysis'], 'transient')

    def test_add_field_current_flow(self):
        w = Agros2DWrapper()

        self.assertRaises(ValueError, w.add_field_current_flow, analysis="eper")
        self.assertRaises(ValueError, w.add_field_current_flow, solver="eper")
        self.assertRaises(ValueError, w.add_field_current_flow, matrix_solver="eper")
        self.assertRaises(ValueError, w.add_field_current_flow, nb_refinements=0.2)
        self.assertRaises(ValueError, w.add_field_current_flow, nb_polyorder=-4)
        self.assertRaises(ValueError, w.add_field_current_flow, adaptivity="eper")

        w.add_field_current_flow(nb_refinement=3)
        self.assertEqual(w.fields['current']['nb_refinement'], 3)


    def test_set_coordinate_type(self):
        w = Agros2DWrapper()

        w.set_coordinate_type('PLANAR')
        self.assertEqual(w.coordinate_type, "planar")

        w.set_coordinate_type('Axisymmetric')
        self.assertEqual(w.coordinate_type, "axisymmetric")

        self.assertRaises(ValueError, w.set_coordinate_type, "eper")


    def test_set_mesh_type(self):
        w = Agros2DWrapper()

        w.set_mesh_type("triangle")
        self.assertEqual(w.mesh_type, "triangle")

        w.set_mesh_type("triangle_quad_fine_division")
        self.assertEqual(w.mesh_type, "triangle_quad_fine_division")

        w.set_mesh_type("triangle_quad_rough_division")
        self.assertEqual(w.mesh_type, "triangle_quad_rough_division")

        w.set_mesh_type("triangle_quad_join")
        self.assertEqual(w.mesh_type, "triangle_quad_join")

        w.set_mesh_type("gmsh_triangle")
        self.assertEqual(w.mesh_type, "gmsh_triangle")

        w.set_mesh_type("gmsh_quad")
        self.assertEqual(w.mesh_type, "gmsh_quad")

        w.set_mesh_type("gmsh_quad_delaunay")
        self.assertEqual(w.mesh_type, "gmsh_quad_delaunay")


        self.assertRaises(ValueError, w.set_mesh_type, "eper")
