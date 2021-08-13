import pygmsh.geo as gmsh
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

    def add_nodes(self):
        # add the adze-nodes to the internal pygmsh-writer object
        # and collects the points into the appropriate format
        for node in self.geometry.nodes:
            self.gmsh_geometry.add_point([node.x, node.y], self.lcar)
            print(self.gmsh_geometry.add_point([node.x, node.y], self.lcar))


    #def add_lines(self):
    #    # add lines to the internal pygmsh object
    #    for line in self.geometry.lines:
    #        start_pt = self.gmsh_geometry.add_point([line.start_pt.x, line.start_pt.y], self.lcar)
    #        end_pt = self.gmsh_geometry.add_point([line.end_pt.x, line.end_pt.y], self.lcar)
    #s        self.gmsh_geometry.add_line(p0=start_pt, p1=end_pt)

    def gmsh_writer(self):
        with gmsh.Geometry() as geom:
            for line in self.geometry.lines:
                start_pt = geom.add_point([line.start_pt.x, line.start_pt.y], self.lcar)
                end_pt = geom.add_point([line.end_pt.x, line.end_pt.y], self.lcar)
                geom.add_line(p0=start_pt, p1=end_pt)
                geom.save_geometry('helo.geo_unrolled')

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
