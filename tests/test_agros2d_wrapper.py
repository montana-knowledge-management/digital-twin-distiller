import os
from unittest import TestCase

from adze_modeler.agros2d_wrapper import Agros2DWrapper


class Agros2DTester(TestCase):
    def test_add_field_electrostatic(self):
        w = Agros2DWrapper()

        self.assertRaises(NotImplementedError, w.add_field, "eper")

        w.add_field("e")
        w.field.set_polynomial_order(3)
        self.assertEqual(w.field.polyorder, 3)

    def test_set_coordinate_type(self):
        w = Agros2DWrapper()

        w.set_coordinate_type("planar")
        self.assertEqual(w.coordinate_type, "planar")

        w.set_coordinate_type("axisymmetric")
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
