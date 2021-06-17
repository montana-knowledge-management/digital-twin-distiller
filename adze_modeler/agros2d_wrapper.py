"""
to execture a script use:
    agros2d_solver -s script.py

execute a script with gui:
    agros2d -s script.py
"""

import sys
from adze_modeler.agros_fields import ElectrostaticField, newline
from adze_modeler.geometry import Geometry


class Agros2DWrapper:
    def __init__(self):
        self.coordinate_type = "planar"
        self.mesh_type = "triangle"

        self.field_e = None  # Electrostatic
        self.field_m = None  # Magnetic
        self.field_h = None  # Heat Flow
        self.field_c = None  # Current Flow

        self.edges = {}

    def add_field(self, type_):

        if type_ == 'e':
            self.field_e = ElectrostaticField()
        else:
            raise NotImplementedError()

    def set_coordinate_type(self, coordinate_type):
        """
        This functions sets the type of the coordinate being used.

        :param coordinate_type: 'planar' or 'axisymmetric'
        """

        if coordinate_type in {"planar", "axisymmetric"}:
            self.coordinate_type = coordinate_type
        else:
            raise ValueError(f'There is no "{coordinate_type}" type of coordinate. ("planar" or "axisymmetric")')

    def set_mesh_type(self, mesh_type):
        """
        This functions sets the type of the mesh being used.

        :param mesh_type: 'triangle', 'triangle_quad_fine_division', 'triangle_quad_rough_division',
        'triangle_quad_join', 'gmsh_triangle', 'gmsh_quad', 'gmsh_quad_delaunay'
        """

        if mesh_type in {
            "triangle",
            "triangle_quad_fine_division",
            "triangle_quad_rough_division",
            "triangle_quad_join",
            "gmsh_triangle",
            "gmsh_quad",
            "gmsh_quad_delaunay",
        }:
            self.mesh_type = mesh_type
        else:
            raise ValueError(f'There is no "{mesh_type}" type of mesh.')

    def add_geometry(self, geo: Geometry):
        for ei in geo.lines:
            key = f'{ei.start_pt.x:.4f}{ei.start_pt.y:.4f}{ei.end_pt.x:.4f}{ei.end_pt.y:.4f}'
            self.edges[key] = (ei.start_pt.x, ei.start_pt.y, ei.end_pt.x, ei.end_pt.y)


    def export(self, out_file):

        sys.stdout = open(out_file, "w")

        print("import agros2d as a2d")
        newline(2)
        print("# problem")
        if self.field_e:
            self.field_e.export()

        if self.field_m:
            self.field_m.export()

        if self.field_h:
            self.field_h.export()

        if self.field_c:
            self.field_c.export()

        sys.stdout.close()
