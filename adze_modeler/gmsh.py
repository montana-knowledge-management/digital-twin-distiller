import pygmsh.geo as gmsh
import adze_modeler.objects as obj
from adze_modeler.geometry import Geometry

"""
The goal of this class is to export the model geometry into a msh file with pygmsh, this mesh file can be 
translated into various formats with the meshio  [1]. 

https://github.com/nschloe/meshio
"""


class GMSHModel:

    def __init__(self, geo, name='dev'):
        self.name = name
        self.boundaries = {}
        self.materials = {}
        self.geometry = geo
        self.metrics = []

        # inner geometry
        self.gmsh_geometry = gmsh.Geometry()

        # sets the
        self.lcar = 5.0  # caracteristic length

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
                                end_point = geom.add_point([edge.start_pt.x, edge.start_pt.y], self.lcar)

                    # in the case of a line
                    if isinstance(edge, obj.Line):
                        line_nr = geom.add_line(p0=start_point, p1=end_point)
                        gmsh_edges.append(line_nr)

                    # circle arcs
                    if isinstance(edge, obj.CircleArc):
                        center_pt = geom.add_point([edge.center_pt.x, edge.center_pt.y], self.lcar)
                        arc_nr = geom.add_circle(start=start_point, center=center_pt, end=end_point)
                        gmsh_edges.append(arc_nr)

                    # bezier curves
                    if isinstance(edge, obj.CubicBezier):
                        center_pt = geom.add_point([edge.center_pt.x, edge.center_pt.y], self.lcar)
                        bezier = geom.add_bspline(start=start_point, center=center_pt, end=end_point)
                        gmsh_edges.append(bezier)


            ll = geom.add_curve_loop(gmsh_edges)
            pl = geom.add_plane_surface(ll)
            geom.save_geometry(file_name + '.geo_unrolled')
            mesh = geom.generate_mesh()
            mesh.write(file_name + ".vtk")

            # for testing the code
            # # plotting out the mesh
            # import pyvista as pv
            # from meshio._helpers import read
            #
            # msh = pv.read(file_name + '.vtk')
            # msh.plot(show_edges=True)
            #
            #
            #
            # print(read(file_name + '.vtk'))
