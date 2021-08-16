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
        self.lcar = 5.0

    def gmsh_writer(self, file_name):
        """
        Writes out the previously defined surfaces from the geo object

        :parameter file_name: the
        """
        with gmsh.Geometry() as geom:
            self.geometry.merge_points()
            surfaces = self.geometry.find_surfaces()

            # the code iterates over the different element types
            for sf in surfaces:
                for edge in sf:
                    # iterate over the line elements
                    if isinstance(edge, obj.Line):
                        # this step can be done during the previous
                        start_pt = geom.add_point([edge.start_pt.x, edge.start_pt.y], self.lcar)
                        end_pt = geom.add_point([edge.end_pt.x, edge.end_pt.y], self.lcar)

                        if edge.id > 0:
                            geom.add_line(p0=start_pt, p1=end_pt)
                        else:
                            geom.add_line(p0=end_pt, p1=start_pt)

                    #
                    if isinstance(edge, obj.CircleArc):
                        # this step can be done during the previous
                        start_pt = geom.add_point([edge.start_pt.x, edge.start_pt.y], self.lcar)
                        center_pt = geom.add_point([edge.center_pt.x, edge.center_pt.y], self.lcar)
                        end_pt = geom.add_point([edge.end_pt.x, edge.end_pt.y], self.lcar)

                        if edge.id > 0:
                            geom.add_circle_arc(start=start_pt, center=center_pt, end=end_pt)
                        else:
                            geom.add_circle(start=end_pt, center=center_pt, end=start_pt)

            geom.save_geometry(file_name + '.geo_unrolled')

# def gmsh_writer(self):
# iterates over all of the closed loops of the geometry

# implemented into the geometry class
# def node_gmsh_point_distance(node, point):
#    dx = node.x - point.x[0]
#    dy = node.y - point.x[1]
#
#    return (dx ** 2.0 + dy ** 2.0) ** 0.5


# def gmsh_writer(nodes, lines, arcs, cubic_beziers):
#     lcar = 5.0
#     epsilon = 1e-6
#     with gmsh.Geometry() as geom:
#         ## add nodes
#         # points = []
#         # for node in nodes:
#         #    temp = geom.add_point([node.x, node.y], lcar)
#         #    points.append(temp)
#
#         # add lines
#         #glines = []
#         #for line in lines:
#         #    for i in range(len(points)):
#         #         if node_gmsh_point_distance(line.start_pt, points[i]) < epsilon:
#         #             start_pt = points[i]
#         #
#         #         if node_gmsh_point_distance(line.end_pt, points[i]) < epsilon:
#         #             end_pt = points[i]
#         #
#         #     temp = geom.add_line(p0=start_pt, p1=end_pt)
#         #     glines.append(temp)
#
#         # add cubic beziers
#         gbeziers = []
#         for cb in cubic_beziers:
#             for i in range(len(points)):
#                 if node_gmsh_point_distance(cb.start_pt, points[i]) < epsilon:
#                     start_pt = points[i]
#                 if node_gmsh_point_distance(cb.end_pt, points[i]) < epsilon:
#                     end_pt = points[i]
#                 if node_gmsh_point_distance(cb.control1, points[i]) < epsilon:
#                     control1 = points[i]
#                 if node_gmsh_point_distance(cb.control2, points[i]) < epsilon:
#                     control2 = points[i]
#
#             temp = geom.add_bspline([start_pt, control1, control2, end_pt])
#             gbeziers.append(temp)
#         # ll = geom.add_curve_loop(glines)
#         # pl = geom.add_plane_surface(ll)
#
#         geom.save_geometry("test.geo_unrolled")
#         # mesh = geom.generate_mesh()
#         # mesh.write("test.vtk")
