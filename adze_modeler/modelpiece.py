from adze_modeler.geometry import Geometry
from adze_modeler.utils import getID


class ModelPiece():

    def __init__(self):
        self.name = None
        self.id = getID()
        self.geom = Geometry()
        self.bbox = [0, 0, 0, 0]
        self.center_point = [0, 0]
        self.label_position = [0, 0] # [0, 1] of bbox  # type: ignore

    def load_piece_from_svg(self, file_name):
        self.geom.import_svg(file_name)
        self._update_bbox()

    def load_piece_from_dxf(self, file_name):
        self.geom.import_dxf(file_name)
        self._update_bbox()

    def translate(self, dx, dy):
        for line_i in self.geom.lines:
            line_i.start_pt.move_xy(dx, dy)
            line_i.end_pt.move_xy(dx, dy)

        self._update_bbox()
        


    def _update_bbox(self):
        minx = min(self.geom.nodes, key=lambda node_i: node_i.x)
        miny = min(self.geom.nodes, key=lambda node_i: node_i.y)
        maxx = max(self.geom.nodes, key=lambda node_i: node_i.x)
        maxy = max(self.geom.nodes, key=lambda node_i: node_i.y)
        self.bbox = [minx, miny, maxx, maxy]
        self.center_point[0] = (maxx-minx) / 2 + minx
        self.center_point[0] = (maxy-miny) / 2 + miny

