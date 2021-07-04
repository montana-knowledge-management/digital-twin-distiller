"""
to execture a script use:
    agros2d_solver -s problem_build_legacy.py

execute a script with gui:
    agros2d -s problem_build_legacy.py
"""
import sys

from adze_modeler.agros_fields import ElectrostaticField
from adze_modeler.agros_fields import HeatFlowField
from adze_modeler.agros_fields import MagneticField
from adze_modeler.agros_fields import newline
from adze_modeler.geometry import Geometry


class Agros2DWrapper:
    def __init__(self):
        self.coordinate_type = "planar"
        self.mesh_type = "triangle"

        # TODO: expand field names
        self.field = None

        self.edges = list()
        self.labels = list()

    def add_field(self, type_):

        if type_ == "e":
            self.field = ElectrostaticField()
        elif type_ == "m":
            self.field = MagneticField()
        elif type_ == "h":
            self.field = HeatFlowField()
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
            self.edges.append((ei.start_pt.x, ei.start_pt.y, ei.end_pt.x, ei.end_pt.y))

    def add_edge(self, start_x, start_y, end_x, end_y, boundary=None):
        self.edges.append((start_x, start_y, end_x, end_y, boundary))

    def add_block_label(self, x, y, name):
        self.labels.append((name, x, y))

    def export(self, out_file):

        # TODO: check if this is portable
        sys.stdout = open(out_file, "w")

        print("import agros2d as a2d")
        newline(2)
        print("# problem")
        print("problem = a2d.problem(clear=True)")
        print(f'problem.coordinate_type = "{self.coordinate_type}"')
        print(f'problem.mesh_type = "{self.mesh_type}"')
        newline(2)
        print("# fields")

        self.field.export()

        newline(2)
        print("# geometry")
        print("geometry = a2d.geometry")
        newline()
        print("# edges")
        for ei in self.edges:
            if ei[4]:
                print(
                    f"geometry.add_edge({ei[0]:.4f}, {ei[1]:.4f}, {ei[2]:.4f}, {ei[3]:.4f}, "
                    f'boundaries={{"{self.field.name}": "{ei[4]}"}})'
                )
            else:
                print(f"geometry.add_edge({ei[0]:.4f}, {ei[1]:.4f}, {ei[2]:.4f}, {ei[3]:.4f})")

        newline()
        print("# block labels")

        for bl_i in self.labels:
            if bl_i[0] in self.field.materials.keys():
                print(
                    f'geometry.add_label({bl_i[1]:.4f}, {bl_i[2]:.4f}, materials = {{"{self.field.name}" : "{bl_i[0]}"}})'
                )
            #
            # elif bl_i[0] in self.field_m.materials.keys():
            #     print(f'geometry.add_label({bl_i[1]:.4f}, {bl_i[2]:.4f}, materials = {{"magnetic" : "{bl_i[0]}"}})')
            #
            # elif bl_i[0] in self.field_c.materials.keys():
            #     # TODO: check is "current" is proper syntax
            #     print(f'geometry.add_label({bl_i[1]:.4f}, {bl_i[2]:.4f}, materials = {{"current" : "{bl_i[0]}"}})')
            #
            # elif bl_i[0] in self.field_h.materials.keys():
            #     print(f'geometry.add_label({bl_i[1]:.4f}, {bl_i[2]:.4f}, materials = {{"heat" : "{bl_i[0]}"}})')

        print("problem.solve()")
        print("a2d.view.zoom_best_fit()")

        sys.stdout.close()
