from copy import copy
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.utils import mirror_point
from adze_modeler.geometry import Geometry
from adze_modeler.metadata import FemmMetadata
from adze_modeler.model import BaseModel
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.objects import Node, Line, CircleArc
import math

def cart2pol(x, y):
    rho = math.hypot(x, y)
    phi = math.atan2(y, x)
    return rho, phi

def pol2cart(rho, phi):
    x = rho * math.cos(math.radians(phi))
    y = rho * math.sin(math.radians(phi))
    return x, y

ORIGIN = Node(0.0, 0.0)
Y = Node(0.0, 1.0)

class BLDCMotor(BaseModel):
    """
    https://www.femm.info/wiki/RotorMotion
    """
    def __init__(self, exportname: str = None):
        super().__init__(exportname)
        self._init_directories()

        # Mesh sizes
        self.msh_size_stator_steel = 1.2
        self.msh_size_rotor_steel = 0.5
        self.msh_size_coils = 1.0
        self.msh_size_air = 1.0
        self.msh_size_airgap = 0.3
        self.msh_size_magnets = 1.0

        # GEOMETRY

        ## ROTOR
        self.r1 = 22.8 / 2   # Rotor Inner Radius
        self.r2 = 50.5 / 2   # Rotor Iron Outer Radius
        self.r3 = 55.1 / 2   # Rotor Outer Radius
        
        ### Magnet
        self.mw = 15.8566    # Magnet Width

        ## AIRGAP
        self.airgap = 0.7
        self.void = 0.3

        ## STATOR
        self.s1 = self.r3 + self.airgap # Stator Inner Radius
        self.s2 = 100 / 2               # Stator Outer Radius

        ## SLOT
        self.w1 = 1.52829
        self.w2 = 3.68306
        self.w3 = 6.8952
        self.w4 = 3.6182

        self.h1 = 0.7
        self.h2 = 0.3707
        self.h3 = 12.1993
        self.h4 = 1.9487


    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False
        femm_metadata.depth = 50

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def add_postprocessing(self):
        pass

    def define_materials(self):
        pass

    def define_boundary_conditions(self):
        pass

    def build_rotor(self):
        g = Geometry()
        alpha1 = math.degrees(math.asin(self.mw * 0.5 / self.r2))
        alpha2 = math.degrees(math.asin(self.mw * 0.5 / self.r3))

        p1l = Node(*pol2cart(self.r1, 90+45/2))
        p2l = Node(*pol2cart(self.r2, 90+45/2))
        p3l = Node(*pol2cart(self.r2, 90+alpha1))
        p4l = Node(*pol2cart(self.r3, 90+alpha2))
        p5l = Node(*pol2cart(self.r3+(self.airgap-self.void) / 2, 90 + 45 / 2))


        p1r = mirror_point(ORIGIN, Y, p1l)
        p2r = mirror_point(ORIGIN, Y, p2l)
        p3r = mirror_point(ORIGIN, Y, p3l)
        p4r = mirror_point(ORIGIN, Y, p4l)
        p5r = mirror_point(ORIGIN, Y, p5l)
        
        g.add_line(Line(p1l, p2l))
        g.add_line(Line(p1r, p2r))
        g.add_line(Line(p3l, p3r))
        g.add_line(Line(p3l, p4l))
        g.add_line(Line(p3r, p4r))
        g.add_line(Line(p2l, p5l))
        g.add_line(Line(p2r, p5r))

        g.add_arc(CircleArc(p1r, ORIGIN, p1l))
        g.add_arc(CircleArc(p3l, ORIGIN, p2l))
        g.add_arc(CircleArc(p2r, ORIGIN, p3r))
        g.add_arc(CircleArc(p4r, ORIGIN, p4l))
        g.add_arc(CircleArc(p5r, ORIGIN, p5l))

        self.geom.merge_geometry(g)

    def build_stator(self):
        g = Geometry()

        q1l = Node(*pol2cart(self.s1 - (self.airgap - self.void) / 2, 90 + 45 / 2))
        q2l = Node(*pol2cart(self.s1, 90 + 45 / 2))
        q3l = Node(*pol2cart(self.s2, 90 + 45 / 2))

        q1r = mirror_point(ORIGIN, Y, q1l)
        q2r = mirror_point(ORIGIN, Y, q2l)
        q3r = mirror_point(ORIGIN, Y, q3l)

        g.add_line(Line(q1l, q2l))
        g.add_line(Line(q2l, q3l))
        g.add_line(Line(q1r, q2r))
        g.add_line(Line(q2r, q3r))

        g.add_arc(CircleArc(q1r, ORIGIN, q1l))
        g.add_arc(CircleArc(q2r, ORIGIN, q2l))
        g.add_arc(CircleArc(q3r, ORIGIN, q3l))

        self.geom.merge_geometry(g)

    def build_slots(self):
        s = ModelPiece("Slot")
        alpha1 = math.degrees(math.asin(self.w1 * 0.5 / self.s1))

        c1l = Node(*pol2cart(self.s1, 90 + alpha1))
        c2l = Node(c1l.x, c1l.y + self.h1)
        c3l = Node(-self.w2 / 2, c2l.y + self.h2)
        c4l = Node(-self.w3 / 2, c3l.y + self.h3)
        c5l = Node(-self.w4 / 2, c4l.y + self.h4)

        c1r = mirror_point(ORIGIN, Y, c1l)
        c2r = mirror_point(ORIGIN, Y, c2l)
        c3r = mirror_point(ORIGIN, Y, c3l)
        c4r = mirror_point(ORIGIN, Y, c4l)
        c5r = mirror_point(ORIGIN, Y, c5l)

        s.geom.add_line(Line(c1l, c2l))
        s.geom.add_line(Line(c2l, c3l))
        s.geom.add_line(Line(c3l, c4l))

        s.geom.add_line(Line(c1r, c2r))
        s.geom.add_line(Line(c2r, c3r))
        s.geom.add_line(Line(c3r, c4r))

        s.geom.add_line(Line(c3l, c3r))

        s.geom.add_arc(CircleArc(c5r, ORIGIN, c5l))
        s.geom.add_arc(CircleArc.from_radius(c5l, c4l, r=1.725))
        s.geom.add_arc(CircleArc.from_radius(c4r, c5r, r=1.725))

        nb_turns = 24
        unitangle = 360 / nb_turns

        s.rotate(alpha=-unitangle)

        for i in range(3):
            si = copy(s)
            si.rotate(alpha=i*unitangle)
            self.geom.merge_geometry(si.geom)


    def build_geometry(self):
        self.build_rotor()
        self.build_stator()
        self.build_slots()
        self.snapshot.add_geometry(self.geom)

if __name__ == "__main__":
    m = BLDCMotor(exportname="dev")
    print(m(devmode=True, cleanup=False))
