from adze_modeler.geometry import Geometry
from adze_modeler.utils import getID


class ModelPiece():

    def __init__(self, name):
        self.name = name
        self.id = getID()
        self.geom = Geometry()
        self.bbox = [0, 0, 0, 0]

    def load_piece_from_svg(self, file_name):
        self.geom.import_svg(str(file_name))
        self._update_bbox()

    def load_piece_from_dxf(self, file_name):
        self.geom.import_dxf(str(file_name))
        self._update_bbox()

    def spawn(self):
        return self.__copy__()

    def translate(self, dx, dy):
        for line_i in self.geom.lines:
            line_i.start_pt.move_xy(dx, dy)
            line_i.end_pt.move_xy(dx, dy)

        self._update_bbox()

    def put(self, x, y):
        deltax = x - self.bbox[0]
        deltay = y - self.bbox[1]
        self.translate(deltax, deltay)

    def _update_bbox(self):
        minx = min(self.geom.nodes, key=lambda node_i: node_i.x).x
        miny = min(self.geom.nodes, key=lambda node_i: node_i.y).y
        maxx = max(self.geom.nodes, key=lambda node_i: node_i.x).x
        maxy = max(self.geom.nodes, key=lambda node_i: node_i.y).y
        self.bbox = [minx, miny, maxx, maxy]

    def __copy__(self):
        piece = ModelPiece(self.name)
        piece.geom.merge_geometry(self.geom)
        piece._update_bbox()
        return piece

