"""
This class realize a layer, where the different elements of the geometry can be stored.
A general geometrical shape can defined by the following objects:
    Nodes (Points), Lines, Circle Arcs, Cubic Bezeirs
"""
import sys
from uuid import uuid4

import svgpathtools

import adze_modeler.objects as obj
import ezdxf
import numpy as np
import svgpathtools as svg
from adze_modeler.utils import getID
from copy import copy


class Geometry:
    def __init__(self):
        self.nodes = []
        self.lines = []
        self.circle_arcs = []
        self.cubic_beziers = []
        self.epsilon = 1.0e-5

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

    def delete_hanging_nodes(self):
        """Delete all nodes, which not part of a another object (Line, Circle, etc)"""
        temp = []
        for node in self.nodes:
            hanging = True
            for line in self.lines:
                if node.id == line.start_pt.id or node.id == line.end_pt.id:
                    hanging = False

            for arc in self.circle_arcs:
                if node.id == arc.start_pt.id or node.id == arc.end_pt.id:
                    hanging = False

            if not hanging:
                temp.append(node)

        del self.nodes
        self.nodes = temp

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

    def merge_lines(self):
        self.merge_points()
        for i in range(len(self.lines) - 1):
            try:
                id1 = self.lines[i].start_pt.id
                id2 = self.lines[i].end_pt.id
            except IndexError:
                pass
            for j in range(len(self.lines) - 1, i, -1):
                id3 = self.lines[j].start_pt.id
                id4 = self.lines[j].end_pt.id
                l = self.lines[j].start_pt.distance_to(self.lines[j].end_pt)

                if l < self.epsilon:
                    del self.lines[j]

                elif {id1, id2} == {id3, id4}:
                    del self.lines[j]

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

        # merge the imported points and coordinates, beca
        self.merge_points()

        return

    def import_svg(self, svg_img, *args):
        """Imports the svg file into a new geo object. The function gives an automatic id to the function

        svg_img: the name of the file, which contains the imported svg image
        return: gives back a new geometry
        """

        # reads the main objects from an svg file
        paths = svg.svg2paths(svg_img)

        # id start from the given number
        id = 0

        for path in paths:
            for seg in path:
                if isinstance(seg, svg.Path):
                    for element in seg:
                        if isinstance(element, svg.Line):
                            start = obj.Node(element.start.real, element.start.imag, id)
                            end = obj.Node(element.end.real, element.end.imag, id + 1)
                            self.add_line(obj.Line(start, end, id + 2))
                            id += 3

                        if isinstance(element, svg.CubicBezier):
                            start = obj.Node(element.start.real, element.start.imag, id)
                            control1 = obj.Node(element.control1.real, element.control1.imag, id + 1)
                            control2 = obj.Node(element.control2.real, element.control2.imag, id + 2)
                            end = obj.Node(element.end.real, element.end.imag, id + 3)
                            self.add_cubic_bezier(obj.CubicBezier(start, control1, control2, end, id + 4))
                            id += 5

        self.merge_points()
        return

    def get_line_intersetions(self, line_1, line_2):
        """
        :param line_1: the first line
        :param line_2: second line

        :returns: None, None or tuple(x, y), None or (x1, y1), (x2, y2)

        https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

        Conditions:
        1. If r × s = 0 and (q − p) × r = 0, then the two lines are collinear.
        2. If r × s = 0 and (q − p) × r ≠ 0, then the two lines are parallel and non-intersecting.
        3. If r × s ≠ 0 and 0 ≤ t ≤ 1 and 0 ≤ u ≤ 1, the two line segments meet at the point p + t r = q + u s.
        4. Otherwise, the two line segments are not parallel but do not intersect.

        """

        x1 = line_1.start_pt.x
        y1 = line_1.start_pt.y
        x2 = line_1.end_pt.x
        y2 = line_1.end_pt.y

        x3 = line_2.start_pt.x
        y3 = line_2.start_pt.y
        x4 = line_2.end_pt.x
        y4 = line_2.end_pt.y

        p1 = None
        p2 = None

        p = np.array([x1, y1])
        r = np.array([x2 - x1, y2 - y1])

        q = np.array([x3, y3])
        s = np.array([x4 - x3, y4 - y3])

        test1 = np.abs(np.cross(r, s))
        test2 = np.abs(np.cross((q - p), r))

        t0 = np.dot(q - p, r) / np.dot(r, r)
        t1 = t0 + np.dot(s, r) / np.dot(r, r)
        t2 = np.dot(p - q, s) / np.dot(s, s)
        t3 = t2 + np.dot(r, s) / np.dot(s, s)

        inrange = lambda x: x > (0 - self.epsilon) and x < (1 + self.epsilon)
        # distance = lambda x, y: np.sqrt((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2)

        if test1 < self.epsilon:

            if test2 < self.epsilon:

                if inrange(t0):
                    p1 = tuple(p + t0 * r)

                if inrange(t1):
                    p2 = tuple(p + t1 * r)

                if inrange(t2):
                    p1 = tuple(q + t2 * s)

                if inrange(t3):
                    p2 = tuple(q + t3 * s)

        else:
            up = (-x1 * y2 + x1 * y3 + x2 * y1 - x2 * y3 - x3 * y1 + x3 * y2) / (
                x1 * y3 - x1 * y4 - x2 * y3 + x2 * y4 - x3 * y1 + x3 * y2 + x4 * y1 - x4 * y2
            )
            tp = (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3) / (
                x1 * y3 - x1 * y4 - x2 * y3 + x2 * y4 - x3 * y1 + x3 * y2 + x4 * y1 - x4 * y2
            )
            if inrange(tp) and inrange(up):
                p1 = tuple(p + tp * r)

        # if (p1 is not None) and (p2 is not None) and (distance(p1, p2) < self.epsilon):
        #     p2 = None

        return p1, p2

    def generate_intersections(self):
        N = len(self.lines)
        newlines = list()
        for i in range(N):
            line_1 = self.lines[i]
            intersections = list()
            distance = lambda a, b: (a.x - b[0]) ** 2 + (a.y - b[1]) ** 2

            for j in range(0, N):  # i+1 is important
                line_2 = self.lines[j]
                p1, p2 = self.get_line_intersetions(line_1, line_2)

                if p1 is not None:
                    if i != j:
                        # plt.scatter(p1[0], p1[1], c="r", marker="o", s=40)
                        intersections.append((distance(line_1.start_pt, p1), *p1))
                        pass

                if p2 is not None:
                    if i != j:
                        # plt.scatter(p2[0], p2[1], c="r", marker="o", s=40)
                        intersections.append((distance(line_1.start_pt, p2), *p2))
                        pass

            intersections.sort(key=lambda ii: ii[0])
            for k in range(len(intersections) - 1):
                start_node = obj.Node(x=intersections[k][1], y=intersections[k][2], id=getID())
                end_node = obj.Node(x=intersections[k + 1][1], y=intersections[k + 1][2], id=getID())
                newlines.append(obj.Line(start_pt=start_node, end_pt=end_node, id=getID()))

        self.nodes.clear()
        self.lines.clear()

        for li in newlines:
            self.add_line(li)

        self.merge_lines()

    def merge_geometry(self, other):

        for li in other.lines:
            otherline = copy(li)
            self.nodes.append(otherline.start_pt)
            self.nodes.append(otherline.end_pt)
            self.lines.append(otherline)

        # for ca in other.circle_arcs:
        #     self.circle_arcs.append(copy(ca))
        #
        # for cb in other.cubic_beziers:
        #     self.cubic_beziers.append(copy(cb))

    def export_geom(self, filename):
        paths = []
        for li in self.lines:
            start_pt = li.start_pt.x + li.start_pt.y * 1j
            end_pt = li.end_pt.x + li.end_pt.y * 1j
            paths.append(svgpathtools.Line(start_pt, end_pt))

        svg.wsvg(paths, filename=str(filename))




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
