"""
This class realize a layer, where the different elements of the geometry can be stored.
A general geometrical shape can defined by the following objects:
    Nodes (Points), Lines, Circle Arcs, Cubic Bezeirs
"""
import sys

import adze_modeler.objects as obj
import ezdxf


class Geometry:
    def __init__(self):
        self.nodes = []
        self.lines = []
        self.circle_arcs = []
        self.cubic_beziers = []
        self.epsilon = 1.0e-5

    def add_object(self):
        return

    def add_node(self, node):
        self.nodes.append(node)

    def add_line(self, line):
        self.lines.append(line)
        # save every start and end points for the geoemtry
        self.nodes.append(line.start_pt)
        self.nodes.append(line.end_pt)

    def add_arc(self, arc):
        self.circle_arcs.append(arc)
        # save every start and end points for the geoemtry
        self.nodes.append(arc.start_pt)
        self.nodes.append(arc.end_pt)
        self.nodes.append(arc.center_pt)

    def add_cubic_bezier(self, cb):
        self.cubic_beziers.append(cb)

        self.nodes.append(cb.start_pt)
        self.nodes.append(cb.control1)
        self.nodes.append(cb.control2)
        self.nodes.append(cb.end_pt)

    def merge_points(self):
        for i in range(len(self.nodes) - 1):
            for j in range(len(self.nodes) - 1, i, -1):
                if self.nodes[i].distance_to(self.nodes[j]) < self.epsilon:

                    # renumber the start/end points of the different shape elements
                    for line in self.lines:
                        if line.start_pt.id == self.nodes[j].id:
                            line.start_pt.id = self.nodes[i].id

                        if line.end_pt.id == self.nodes[j].id:
                            line.end_pt.id = self.nodes[i].id

                    del self.nodes[j]
            # if node_1.distance_to(node_2) < self.epsilon:
            #    print(i)

    def meshi_it(self, mesh_strategy):
        mesh = mesh_strategy(self.nodes, self.lines, self.circle_arcs, self.cubic_beziers)
        return mesh

    def __repr__(self):
        msg = ""
        msg += "\n Nodes:       \n -----------------------\n"
        for node in self.nodes:
            msg += str(node) + "\n"

        msg += "\n Lines:       \n -----------------------\n"
        for line in self.lines:
            msg += str(line) + "\n"

        msg += "\n Circle Arcs: \n -----------------------\n"
        for arc in self.circle_arcs:
            msg += str(arc) + " \n"

        msg += "\n CubicBezier: \n -----------------------\n"
        for cubicbezier in self.cubic_beziers:
            msg += str(cubicbezier) + "\n"

        return msg

    def import_dxf(self, dxf_file):
        try:
            doc = ezdxf.readfile(dxf_file)
        except OSError:
            print("Not a DXF file or a generic I/O error.")
            sys.exit(1)
        except ezdxf.DXFStructureError:
            print("Invalid or corrupted DXF file.")
            sys.exit(2)

        # iterate over all entities in modelspace
        # id start from the given number
        id = 0

        msp = doc.modelspace()
        for e in msp:
            if e.dxftype() == "LINE":
                start = obj.Node(e.dxf.start[0], e.dxf.start[1], id)
                end = obj.Node(e.dxf.end[0], e.dxf.end[1], id + 1)
                self.add_line(obj.Line(start, end, id + 2))
                id += 3

            if e.dxftype() == "ARC":
                start = obj.Node(e.start_point.x, e.start_point.y, id)
                end = obj.Node(e.end_point.x, e.end_point.y, id + 1)
                center = obj.Node(e.dxf.center[0], e.dxf.center[0], id + 2)
                self.add_arc(obj.CircleArc(start, center, end, id + 3))
                id += 4

            if e.dxftype() == "POLYLINE":
                print(e.__dict__)
        return


# def node_gmsh_point_distance(node, point):
#     dx = node.x - point.x[0]
#     dy = node.y - point.x[1]
#
#     return (dx ** 2. + dy ** 2.) ** 0.5


# def gmsh_writer(nodes, lines, arcs, cubic_beziers):
#     lcar = 5.
#     epsilon = 1e-6
#     with gmsh.Geometry() as geom:
#         # add nodes
#         points = []
#         for node in nodes:
#             temp = geom.add_point([node.x, node.y], lcar)
#             # temp._id = node.id
#             points.append(temp)
#
#         # add lines
#         glines = []
#         for line in lines:
#             for i in range(len(points)):
#                 if node_gmsh_point_distance(line.start_pt, points[i]) < epsilon:
#                     start_pt = points[i]
#
#                 if node_gmsh_point_distance(line.end_pt, points[i]) < epsilon:
#                     end_pt = points[i]
#
#             temp = geom.add_line(p0=start_pt, p1=end_pt)
#             glines.append(temp)
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
