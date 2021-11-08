"""
This class realize a layer, where the different elements of the geometry can be stored.

A general geometrical shape can defined by the following objects:
    Nodes (Points), Lines, Circle Arcs, Cubic Bezeirs
"""
import os
import sys
from copy import copy, deepcopy

import ezdxf
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import svgpathtools as svg

import digital_twin_distiller.objects as obj
from digital_twin_distiller.utils import getID


class Geometry:
    def __init__(self):
        self.nodes = []
        self.lines = []
        self.circle_arcs = []
        self.cubic_beziers = []
        self.epsilon = 1.0e-5

    def add_node(self, node):
        # self.nodes.append(copy(node))
        self.append_node(node)
        return node

    def add_line(self, line):
        # save every start and end points for the geoemtry
        line.start_pt = self.append_node(line.start_pt)
        line.end_pt = self.append_node(line.end_pt)

        self.lines.append(line)

    def add_arc(self, arc):
        # save every start and end points for the geoemtry if they are not exists
        arc.start_pt = self.append_node(arc.start_pt)
        arc.end_pt = self.append_node(arc.end_pt)

        self.circle_arcs.append(arc)

    def add_cubic_bezier(self, cb):
        # save every start and end points for the geoemtry if they are not exists
        cb.start_pt = self.append_node(cb.start_pt)
        cb.end_pt = self.append_node(cb.end_pt)
        self.cubic_beziers.append(cb)

    def add_rectangle(self, r: obj.Rectangle):
        p = list(r)
        self.add_line(obj.Line(obj.Node(p[0], p[1]), obj.Node(p[2], p[3])))
        self.add_line(obj.Line(obj.Node(p[2], p[3]), obj.Node(p[4], p[5])))
        self.add_line(obj.Line(obj.Node(p[4], p[5]), obj.Node(p[6], p[7])))
        self.add_line(obj.Line(obj.Node(p[6], p[7]), obj.Node(p[0], p[1])))

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

    def append_node(self, new_node):
        """Appends the node to the node list only if its not exists, gives back that node object """
        for i in range(len(self.nodes)):
            if self.nodes[i].distance_to(new_node) < self.epsilon:
                return self.nodes[i]

        self.nodes.append(new_node)
        return new_node

    def merge_points(self):
        for i in range(len(self.nodes) - 1):
            for j in range(len(self.nodes) - 1, i, -1):
                if self.nodes[i].distance_to(self.nodes[j]) < self.epsilon:

                    # renumber the start/end points of the different shape elements
                    # 1/ lines
                    for line in self.lines:
                        if line.start_pt.id == self.nodes[j].id:
                            line.start_pt.id = self.nodes[i].id

                        if line.end_pt.id == self.nodes[j].id:
                            line.end_pt.id = self.nodes[i].id

                    # 2/ arcs
                    for arcs in self.circle_arcs:
                        if arcs.start_pt.id == self.nodes[j].id:
                            arcs.start_pt.id = self.nodes[i].id

                        if arcs.end_pt.id == self.nodes[j].id:
                            arcs.end_pt.id = self.nodes[i].id

                    # 3/ bezier curves
                    for cb in self.cubic_beziers:
                        if cb.start_pt.id == self.nodes[j].id:
                            cb.start_pt.id = self.nodes[i].id

                        if cb.control1.id == self.nodes[j].id:
                            cb.control1.id = self.nodes[i].id

                        if cb.control2.id == self.nodes[j].id:
                            cb.control2.id = self.nodes[i].id

                        if cb.end_pt.id == self.nodes[j].id:
                            cb.end_pt.id = self.nodes[i].id

                    del self.nodes[j]

    # Todo: bezierre és circle arcra be kellene fejezni
    def merge_lines(self):
        lines = self.lines.copy()
        self.lines.clear()

        for li in lines:
            if li not in self.lines:
                self.add_line(li)

        self.merge_points()

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
            doc = ezdxf.readfile(str(dxf_file))
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
                start = obj.Node(e.dxf.start[0], e.dxf.start[1])
                end = obj.Node(e.dxf.end[0], e.dxf.end[1])
                self.add_line(obj.Line(start, end))
                id += 3

            if e.dxftype() == "ARC":
                start = obj.Node(e.start_point.x, e.start_point.y)
                end = obj.Node(e.end_point.x, e.end_point.y)
                center = obj.Node(e.dxf.center[0], e.dxf.center[1])

                self.add_arc(obj.CircleArc(start, center, end))

                id += 4

            if e.dxftype() == "POLYLINE":
                print(e.__dict__)

        # merge the imported points and coordinates, beca
        self.merge_points()

        return

    @staticmethod
    def casteljau(bezier: obj.CubicBezier):
        """
        Gets a Bezier object and makes only one Casteljau's iteration step on it without the recursion.

        The algorithm splits the bezier into two, smaller parts denoted by r is the 'right-sides' and l denotes the
        'left sided' one. The algorithm is limited to use cubic beziers only.

        :return: 2 bezier objects, the right and the left one

        """
        # calculating the mid point [m]
        m = (bezier.control1 + bezier.control2) * 0.5

        l0 = bezier.start_pt
        r3 = bezier.end_pt

        l1 = (bezier.start_pt + bezier.control1) * 0.5
        r2 = (bezier.control2 + bezier.end_pt) * 0.5

        l2 = (l1 + m) * 0.5
        r1 = (r2 + m) * 0.5

        l3 = (l2 + r1) * 0.5
        r0 = l3

        r = obj.CubicBezier(start_pt=r0, control1=r1, control2=r2, end_pt=r3)
        l = obj.CubicBezier(start_pt=l0, control1=l1, control2=l2, end_pt=l3)

        return r, l

    def export_svg(self, file_name="output.svg"):
        """
        Creates an svg image from the geometry objects.
        """

        # every object handled as a separate path

        obj.Line
        paths = []
        colors = []
        # exports the lines
        for seg in self.lines:
            path = svg.Path()
            path.append(
                svg.Line(
                    complex(seg.start_pt.x, seg.start_pt.y).conjugate(),
                    complex(seg.end_pt.x, seg.end_pt.y).conjugate(),
                )
            )
            paths.append(path)
            colors.append("blue")

        # export the circle arcs
        # TODO: be kellene fejezni
        for arc in self.circle_arcs:
            path = svg.Path()
            path.append(
                svg.Line(
                    complex(arc.start_pt.x, arc.start_pt.y).conjugate(),
                    complex(arc.end_pt.x, arc.end_pt.y).conjugate(),
                )
            )
            paths.append(path)
            colors.append("blue")

        for cb in self.cubic_beziers:
            path = svg.Path()
            p1 = complex(cb.start_pt.x, cb.start_pt.y).conjugate()
            p2 = complex(cb.end_pt.x, cb.end_pt.y).conjugate()
            c1 = complex(cb.control1.x, cb.control1.y).conjugate()
            c2 = complex(cb.control2.x, cb.control2.y).conjugate()
            path.append(svg.CubicBezier(p1, c1, c2, p2))
            paths.append(path)
            colors.append("blue")

        work_dir = os.getcwd()
        svg.wsvg(paths, colors=colors, svgwrite_debug=True, filename=work_dir + "/" + file_name)

    def import_svg(self, svg_img, *args):
        """Imports the svg file into a new geo object. The function gives an automatic id to the function

        svg_img: the name of the file, which contains the imported svg image
        return: gives back a new geometry
        """

        # reads the main objects from an svg file
        paths = svg.svg2paths(str(svg_img))

        # id start from the given number
        id = 0

        for path in paths:
            for seg in path:
                if isinstance(seg, svg.Path):
                    for element in seg:
                        if isinstance(element, svg.Line):
                            p1 = element.start.conjugate()
                            p2 = element.end.conjugate()
                            start = obj.Node(p1.real, p1.imag)
                            end = obj.Node(p2.real, p2.imag)
                            self.add_line(obj.Line(start, end))
                            id += 3

                        if isinstance(element, svg.CubicBezier):
                            s1 = element.start.conjugate()
                            s2 = element.end.conjugate()
                            c1 = element.control1.conjugate()
                            c2 = element.control2.conjugate()

                            start = obj.Node(s1.real, s1.imag, id)
                            control1 = obj.Node(c1.real, c1.imag, id + 1)
                            control2 = obj.Node(c2.real, c2.imag, id + 2)
                            end = obj.Node(s2.real, s2.imag, id + 3)
                            self.add_cubic_bezier(obj.CubicBezier(start, control1, control2, end, id + 4))
                            id += 5

                        if isinstance(element, svg.Arc):
                            # pass
                            p1 = element.start.conjugate()
                            p2 = element.center.conjugate()
                            p3 = element.end.conjugate()
                            start = obj.Node(p1.real, p1.imag)
                            center = obj.Node(p2.real, p2.imag)
                            end = obj.Node(p3.real, p3.imag)
                            self.add_arc(obj.CircleArc(start, center, end))

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
        """Todo: generate intersections"""
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

                if p2 is not None:
                    if i != j:
                        intersections.append((distance(line_1.start_pt, p2), *p2))

            intersections.sort(key=lambda ii: ii[0])
            for k in range(len(intersections) - 1):
                start_node = obj.Node(x=intersections[k][1], y=intersections[k][2], id=getID())
                end_node = obj.Node(
                    x=intersections[k + 1][1],
                    y=intersections[k + 1][2],
                    id=getID(),
                )
                newlines.append(obj.Line(start_pt=start_node, end_pt=end_node, id=getID()))

        self.nodes.clear()
        self.lines.clear()

        for li in newlines:
            if li.start_pt.distance_to(li.end_pt) > self.epsilon:
                self.add_line(li)

        self.merge_lines()

    def merge_geometry(self, other):

        for ni in other.nodes:
            self.add_node(copy(ni))

        for li in other.lines:
            self.add_line(copy(li))

        for i, ca in enumerate(other.circle_arcs):
            self.add_arc(copy(ca))

        for cb in other.cubic_beziers:
            self.add_cubic_bezier(copy(cb))

    # Todo: szerintem ez az export svg reinkarnációja lehet
    # def export_geom(self, filename):
    #     paths = []
    #     for li in self.lines:
    #         start_pt = li.start_pt.x + li.start_pt.y * 1j
    #         end_pt = li.end_pt.x + li.end_pt.y * 1j
    #         paths.append(svgpathtools.Line(start_pt, end_pt))
    #
    #     for bz in self.cubic_beziers:
    #         start_pt = complex(*bz.start_pt)
    #         control1 = complex(*bz.control1)
    #         control2 = complex(*bz.control2)
    #         end_pt = complex(*bz.end_pt)
    #         paths.append(svgpathtools.CubicBezier(start_pt, control1, control2, end_pt))
    #
    #     svg.wsvg(paths, filename=str(filename))

    def __copy__(self):
        g = Geometry()
        g.merge_geometry(self)
        return g

    def find_edges(self, nodes: list):
        """Search for the edges with the given direction"""

        surface = []
        # we are looking for a closed loop, therefore the first and the last item should create an edge
        nodes.append(nodes[0])
        for i, node in enumerate(nodes[:-1]):
            # lines
            for line in self.lines:
                if line.end_pt.id == node and line.start_pt.id == nodes[i + 1]:
                    # direction: from end -> to start
                    temp = deepcopy(line)
                    temp.id = -temp.id
                    surface.append(temp)

                if line.start_pt.id == node and line.end_pt.id == nodes[i + 1]:
                    # direction: from start -> to end
                    surface.append(deepcopy(line))

            # arcs
            for arc in self.circle_arcs:
                if arc.end_pt.id == node and arc.start_pt.id == nodes[i + 1]:
                    # direction: from end -> to start
                    temp = deepcopy(arc)
                    temp.id = -temp.id
                    surface.append(temp)

                if arc.start_pt.id == node and arc.end_pt.id == nodes[i + 1]:
                    # direction: from start -> to end
                    surface.append(deepcopy(arc))

            # cubic bezier
            for cb in self.cubic_beziers:
                if cb.end_pt.id == node and cb.start_pt.id == nodes[i + 1]:
                    # direction: from end -> to start
                    temp = deepcopy(cb)
                    temp.id = -temp.id
                    surface.append(temp)

                if cb.start_pt.id == node and cb.end_pt.id == nodes[i + 1]:
                    # direction: from start -> to end
                    surface.append(deepcopy(cb))

        return surface

    def find_surfaces(self):
        """Builds a networkx graph from the given geometry and search for the closed surfaces with networkx."""

        Graph = nx.Graph()

        # add edges to the graph from the different entities: lines, circles, cubic_bezier
        for line in self.lines:
            Graph.add_edge(line.start_pt.id, line.end_pt.id)

        for circles in self.circle_arcs:
            Graph.add_edge(circles.start_pt.id, circles.end_pt.id)

        for cb in self.cubic_beziers:
            Graph.add_edge(cb.start_pt.id, cb.end_pt.id)

        closed_loops = nx.cycle_basis(Graph)

        surfaces = []

        for loop in closed_loops:
            surface = self.find_edges(loop)
            surfaces.append(surface)

        return surfaces

    def plot_connection_graph(self, debug=False):
        """Plots the connection graph of the given task."""
        Graph = nx.Graph()

        # add edges to the graph from the different entities: lines, circles, cubic_bezier
        for line in self.lines:
            Graph.add_edge(line.start_pt.id, line.end_pt.id)

        for circles in self.circle_arcs:
            Graph.add_edge(circles.start_pt.id, circles.end_pt.id)

        for cb in self.cubic_beziers:
            Graph.add_edge(cb.start_pt.id, cb.end_pt.id)

        nx.draw_networkx(Graph, with_labels=False)

        if not debug:
            # Set margins for the axes so that nodes aren't clipped
            ax = plt.gca()
            ax.margins(0.20)
            plt.axis("off")
            plt.show()
