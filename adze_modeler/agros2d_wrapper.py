"""
to execture a script use:
    agros2d_solver -s script.py

execute a script with gui:
    agros2d -s script.py
"""

import sys


class Agros2DWrapper:
    def __init__(self):
        self.coordinate_type = "planar"
        self.mesh_type = "triangle"

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

    def add_boundary_condition(self, name, value, bctype, field):
        """
        Define a new boundary condition.

        :param name: name of the boundary condition
        :param value: The value of the boundary condition
        :param bctype: 'dirichlet', 'neumann'
        :param field: 'electrostatic', 'magnetic', 'heat', 'current'
        """
        if name not in self.boundary_conditions.keys():
            self.boundary_conditions[name] = (field, bctype, value)
        else:
            raise ValueError(f'There is a boundary condition called "{name}".')

    def add_material(self, name, **kwargs):
        pass


    def export(self, out_file):
        # TODO: Add field keywords correcting function

        sys.stdout = open(out_file, 'w')
        newline = lambda n: print('\n'*(n-1))

        print("import agros2d as a2d")
        newline(2)
        print("# problem")

        sys.stdout.close()