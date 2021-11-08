import gmsh as std_gmsh
import meshio
import pygmsh.geo as gmsh

import digital_twin_distiller.objects as obj
from digital_twin_distiller.geometry import Geometry

"""
The goal of this class is to export the model geometry into a msh file with pygmsh, this mesh file can be
translated into various formats with the meshio  [1].

https://github.com/nschloe/meshio

useful documentation for the usage of gmsh - pygmsh codes:
http://jsdokken.com/converted_files/tutorial_pygmsh.html
"""


class GMSHModel:
    def __init__(self, geo, name="dev", msh_format=".msh"):
        self.name = name
        self.boundaries = {}  # this should be defined
        self.boundary_queue_gmsh = {}  # gmsh renumbers the different nodes and
        self.label_queue = []
        self.materials = {}
        self.geometry = geo
        self.metrics = []

        # inner geometry
        self.gmsh_geometry = gmsh.Geometry()

        # sets the
        self.lcar = 5  # characteristic length
        self.msh_format = msh_format
        self.dim = 2  # dimension of the mesh

    @staticmethod
    def export_mesh_to_dolphin(mesh, cell_type, prune_z=False):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        out_mesh = meshio.Mesh(
            points=mesh.points,
            cells={cell_type: cells},
            cell_data={"name_to_read": [cell_data]},
        )
        if prune_z:
            out_mesh.prune_z_0()
        return out_mesh

    def gmsh_writer(self, file_name):
        """
        Writes out the previously defined surfaces from the geo object

        :parameter file_name: the
        """
        gmsh_edges = []  # the id numbers for the gmsh edges

        with gmsh.Geometry() as geom:
            self.geometry.merge_points()
            surfaces = self.geometry.find_surfaces()

            # the code iterates over the different element types
            for sf in surfaces:

                # firstly, we have to build a closed loop from the edges of the surface, this closed loop should be
                # a directed graph, therefore it is important to write out the lines in the right order
                closed_loop = []
                start_point = None
                end_point = None
                for index, edge in enumerate(sf):

                    # firstly, the code ordering the lines into the right order, to form a directed closed loop
                    if not start_point:
                        if edge.id > 0:
                            start_point = geom.add_point([edge.start_pt.x, edge.start_pt.y], self.lcar)
                            end_point = geom.add_point([edge.end_pt.x, edge.end_pt.y], self.lcar)

                        else:
                            start_point = geom.add_point([edge.end_pt.x, edge.end_pt.y], self.lcar)
                            end_point = geom.add_point([edge.start_pt.x, edge.start_pt.y], self.lcar)

                        first_point = start_point
                    else:
                        start_point = end_point
                        # closing the loop if this is the final edge in the list
                        if index == len(sf) - 1:
                            end_point = first_point
                        else:
                            if edge.id > 0:
                                end_point = geom.add_point([edge.end_pt.x, edge.end_pt.y], self.lcar)
                            else:
                                end_point = geom.add_point(
                                    [edge.start_pt.x, edge.start_pt.y],
                                    self.lcar,
                                )

                    # in the case of a line
                    if isinstance(edge, obj.Line):
                        line_nr = geom.add_line(p0=start_point, p1=end_point)
                        gmsh_edges.append(line_nr)

                    # circle arcs
                    if isinstance(edge, obj.CircleArc):
                        center_pt = geom.add_point([edge.center_pt.x, edge.center_pt.y], self.lcar)
                        # arc_nr = geom.add_circle(start=start_point, center=center_pt, end=end_point)
                        arc_nr = geom.add_circle_arc(start=start_point, center=center_pt, end=end_point)
                        gmsh_edges.append(arc_nr)

                    # bezier curves
                    if isinstance(edge, obj.CubicBezier):
                        control1 = geom.add_point([edge.control1.x, edge.control1.y], self.lcar)
                        control2 = geom.add_point([edge.control2.x, edge.control2.y], self.lcar)
                        bezier = geom.add_bspline(
                            control_points=[
                                start_point,
                                control1,
                                control2,
                                end_point,
                            ]
                        )
                        gmsh_edges.append(bezier)

                    # the number of the boundaries should be renumbered to get the
                    for key, val in self.boundaries.items():
                        if abs(edge.id) in val:
                            if key in self.boundary_queue_gmsh:
                                self.boundary_queue_gmsh[key].append(gmsh_edges[-1])
                            else:
                                self.boundary_queue_gmsh[key] = [gmsh_edges[-1]]

            ll = geom.add_curve_loop(gmsh_edges)
            pl = geom.add_plane_surface(ll)

            # physical surfaces define the name of the applied materials
            geom.add_physical(pl, label="material")

            # define the boundary condition for the latest edge
            for key, val in self.boundary_queue_gmsh.items():
                # gmsh define the boundaries in the transposed order
                geom.add_physical(val, label=key)

            # geom.save_geometry(file_name + '.geo_unrolled') physical domain saving not working with it
            mesh = geom.generate_mesh(dim=self.dim)
            std_gmsh.write(file_name + ".msh")
            std_gmsh.write(file_name + ".geo_unrolled")
            std_gmsh.clear()

            # # create a mesh for dolphin
            # mesh_from_file = meshio.read(file_name + ".msh")
            # # line_mesh = self.create_mesh(mesh_from_file, "line", prune_z=True)
            # # meshio.write(file_name + "_bounds_" + ".xdmf", line_mesh)
            #
            # cell_mesh = self.create_mesh(mesh_from_file, "triangle", prune_z=True)
            # meshio.write(file_name + ".xdmf", cell_mesh)
            #
            # # if self.msh_format != '.msh':
            # #    mesh_from_file = meshio.read(file_name + ".msh")
            # #    meshio.write(self.msh_format, mesh_from_file)
            # # mesh.write(file_name + ".xdmf")  # + ".vtk")
