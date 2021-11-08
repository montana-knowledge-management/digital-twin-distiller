import operator
from math import pi

from digital_twin_distiller.geometry import Geometry
from digital_twin_distiller.objects import Node
from digital_twin_distiller.utils import getID, get_short_id, mirror_point


class ModelPiece:
    def __init__(self, name):
        self.name = name
        self.id = getID()
        self.geom = Geometry()
        self.bbox = [0, 0, 0, 0]

        # extreme points
        self.left = None
        self.right = None
        self.upper = None
        self.lower = None

    def load_piece_from_svg(self, file_name):
        self.geom.import_svg(str(file_name))
        self.update_bbox()

    def load_piece_from_dxf(self, file_name):
        self.geom.import_dxf(str(file_name))
        self.update_bbox()

    def spawn(self):
        return self.__copy__()

    def translate(self, dx, dy):
        updated = set()
        for i, node_i in enumerate(self.geom.nodes):
            self.geom.nodes[i].move_xy(dx, dy)
            updated.add(get_short_id(node_i))

        for i, li in enumerate(self.geom.lines):
            id_start = get_short_id(li.start_pt)
            id_end = get_short_id(li.end_pt)
            if id_start not in updated:
                self.geom.lines[i].start_pt.move_xy(dx, dy)
                updated.add(id_start)

            if id_end not in updated:
                self.geom.lines[i].end_pt.move_xy(dx, dy)
                updated.add(id_end)

        for i, ai in enumerate(self.geom.circle_arcs):
            id_start = get_short_id(ai.start_pt)
            id_center = get_short_id(ai.center_pt)
            id_apex = get_short_id(ai.apex_pt)
            id_end = get_short_id(ai.end_pt)

            if id_start not in updated:
                self.geom.circle_arcs[i].start_pt.move_xy(dx, dy)
                updated.add(id_start)

            if id_center not in updated:
                self.geom.circle_arcs[i].center_pt.move_xy(dx, dy)
                updated.add(id_center)

            if id_apex not in updated:
                self.geom.circle_arcs[i].apex_pt.move_xy(dx, dy)
                updated.add(id_apex)

            if id_end not in updated:
                self.geom.circle_arcs[i].end_pt.move_xy(dx, dy)
                updated.add(id_end)

        self.update_bbox()

    def put(self, x, y, bbox_ref="lower-left"):
        if bbox_ref == "lower-left":
            deltax = x - self.bbox[0]
            deltay = y - self.bbox[1]

        elif bbox_ref == "lower-right":
            deltax = x - self.bbox[2]
            deltay = y - self.bbox[1]

        elif bbox_ref == "upper-left":
            deltax = x - self.bbox[0]
            deltay = y - self.bbox[3]

        elif bbox_ref == "upper-right":
            deltax = x - self.bbox[2]
            deltay = y - self.bbox[3]

        elif bbox_ref == "upper":
            deltax = x - self.upper.x
            deltay = y - self.upper.y

        elif bbox_ref == "lower":
            deltax = x - self.lower.x
            deltay = y - self.lower.y

        elif bbox_ref == "right":
            deltax = x - self.right.x
            deltay = y - self.right.y

        elif bbox_ref == "left":
            deltax = x - self.left.x
            deltay = y - self.left.y
        else:
            raise ValueError()

        self.translate(deltax, deltay)

    def mirror(self, p1=(0, 0), p2=(0, 1)):
        """
        This function mirrors all the geometry point on a line defined by p1 and p2.
        """
        p1 = Node(*p1)
        p2 = Node(*p2)
        for i in range(len(self.geom.nodes)):
            self.geom.nodes[i] = mirror_point(p1, p2, self.geom.nodes[i])

        for i in range(len(self.geom.circle_arcs)):
            # when mirroring the start and end points are going to be swapped to preserve the arc direction
            start = mirror_point(p1, p2, self.geom.circle_arcs[i].start_pt)
            end = mirror_point(p1, p2, self.geom.circle_arcs[i].end_pt)
            self.geom.circle_arcs[i].start_pt = end
            self.geom.circle_arcs[i].center_pt = mirror_point(p1, p2, self.geom.circle_arcs[i].center_pt)
            self.geom.circle_arcs[i].apex_pt = mirror_point(p1, p2, self.geom.circle_arcs[i].apex_pt)
            self.geom.circle_arcs[i].end_pt = start

        for i in range(len(self.geom.lines)):
            self.geom.lines[i].start_pt = mirror_point(p1, p2, self.geom.lines[i].start_pt)
            self.geom.lines[i].end_pt = mirror_point(p1, p2, self.geom.lines[i].end_pt)

    def rotate(self, ref_point=(0, 0), alpha=0):
        """
        Rotate all points of the modelpiece around the reference point with alpha degrees.
        """

        ref_point = Node(*ref_point)
        alpha = alpha * pi / 180
        for i in range(len(self.geom.nodes)):
            self.geom.nodes[i] = self.geom.nodes[i].rotate_about(ref_point, alpha)

        for i in range(len(self.geom.circle_arcs)):
            self.geom.circle_arcs[i].start_pt = self.geom.circle_arcs[i].start_pt.rotate_about(ref_point, alpha)
            self.geom.circle_arcs[i].center_pt = self.geom.circle_arcs[i].center_pt.rotate_about(ref_point, alpha)
            self.geom.circle_arcs[i].apex_pt = self.geom.circle_arcs[i].apex_pt.rotate_about(ref_point, alpha)
            self.geom.circle_arcs[i].end_pt = self.geom.circle_arcs[i].end_pt.rotate_about(ref_point, alpha)

        for i in range(len(self.geom.lines)):
            self.geom.lines[i].start_pt = self.geom.lines[i].start_pt.rotate_about(ref_point, alpha)
            self.geom.lines[i].end_pt = self.geom.lines[i].end_pt.rotate_about(ref_point, alpha)

        self.update_bbox()

    def scale(self, sx, sy):
        updated = set()
        for i in range(len(self.geom.nodes)):
            self.geom.nodes[i].x *= sx
            self.geom.nodes[i].y *= sy
            updated.add(get_short_id(self.geom.nodes[i]))

        for i in range(len(self.geom.lines)):
            id_start = get_short_id(self.geom.lines[i].start_pt)
            id_end = get_short_id(self.geom.lines[i].end_pt)

            if id_start not in updated:
                self.geom.lines[i].start_pt.x *= sx
                self.geom.lines[i].start_pt.y *= sy
                updated.add(id_start)

            if id_end not in updated:
                self.geom.lines[i].end_pt.x *= sx
                self.geom.lines[i].end_pt.y *= sy
                updated.add(id_end)

        for i in range(len(self.geom.circle_arcs)):
            id_start = get_short_id(self.geom.circle_arcs[i].start_pt)
            id_center = get_short_id(self.geom.circle_arcs[i].center_pt)
            id_apex = get_short_id(self.geom.circle_arcs[i].apex_pt)
            id_end = get_short_id(self.geom.circle_arcs[i].end_pt)

            if id_start not in updated:
                self.geom.circle_arcs[i].start_pt.x *= sx
                self.geom.circle_arcs[i].start_pt.y *= sy
                updated.add(id_start)

            if id_center not in updated:
                self.geom.circle_arcs[i].center_pt.x *= sx
                self.geom.circle_arcs[i].center_pt.y *= sy
                updated.add(id_center)

            if id_apex not in updated:
                self.geom.circle_arcs[i].apex_pt.x *= sx
                self.geom.circle_arcs[i].apex_pt.y *= sy
                updated.add(id_apex)

            if id_end not in updated:
                self.geom.circle_arcs[i].end_pt.x *= sx
                self.geom.circle_arcs[i].end_pt.y *= sy
                updated.add(id_end)

    def update_bbox(self):
        minx = min(self.geom.nodes, key=lambda node_i: node_i.x).x
        miny = min(self.geom.nodes, key=lambda node_i: node_i.y).y
        maxx = max(self.geom.nodes, key=lambda node_i: node_i.x).x
        maxy = max(self.geom.nodes, key=lambda node_i: node_i.y).y
        self.bbox = [minx, miny, maxx, maxy]
        self._update_extreme_points()

    def _update_extreme_points(self):
        self.upper = max(self.geom.nodes, key=operator.attrgetter("y"))
        self.lower = min(self.geom.nodes, key=operator.attrgetter("y"))
        self.right = max(self.geom.nodes, key=operator.attrgetter("x"))
        self.left = min(self.geom.nodes, key=operator.attrgetter("x"))

    def __copy__(self):
        piece = ModelPiece(self.name)
        piece.geom.merge_geometry(self.geom)
        piece.update_bbox()
        return piece
