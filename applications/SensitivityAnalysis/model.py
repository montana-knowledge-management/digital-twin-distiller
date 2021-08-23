import math
from copy import copy
from pathlib import Path
from shutil import rmtree

import numpy as np
from math import pi, sin, cos, tan

from adze_modeler.snapshot import Snapshot
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.model import BaseModel
from adze_modeler.metadata import FemmMetadata
from adze_modeler.platforms.femm import Femm
from adze_modeler.material import Material
from adze_modeler.boundaries import AntiPeriodicBoundaryCondition, DirichletBoundaryCondition, AntiPeriodicAirGap
from adze_modeler.objects import Line, Node, CircleArc, Rectangle
import numpy as np
import scipy


class PriusMotor(BaseModel):
    """docstring for PriusMotor"""

    def __init__(self, rotorangle=0, exportname=None):
        super(PriusMotor, self).__init__(exportname=exportname)

        self.rotorangle = rotorangle

        self.S1 = 134.62
        self.S2 = 80.95
        self.S4 = 1.92953
        self.R1 = 80.2
        self.R3 = 140
        self.R4 = 6.5
        self.R5 = 18.9
        self.R6 = 70.08
        self.R7 = 55.5

        self.airgap = 0.5

        self.label_queue = []
        self.boundary_queue = []
        self.boundary_arc_queue = []

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False
        femm_metadata.depth = 84
        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def add_postprocessing(self):
        entities = [
            (0, self.R6 - 0.1),
            (0, self.R1 - 0.1),
            (-10, 70),
            (10, 70),
            (-20, 75),
            (20, 75)
        ]
        self.snapshot.add_postprocessing('integration', entities, "Torque")

    def define_materials(self):
        air = Material("air")
        core = Material("core")
        magnet = Material("magnet")

        awg2 = Material("20 AWG")
        awg2.lamination_type = "magnetwire"
        awg2.diameter = 0.812049969500513
        awg2.conductivity = 58e6

        coil_area = 0.000143002
        J = 25 / coil_area
        phase_A = copy(awg2)
        phase_A.name = "A+"
        phase_A.Je = - J / 2

        phase_B = copy(awg2)
        phase_B.name = "B-"
        phase_B.Je = J

        phase_C = copy(awg2)
        phase_C.name = "C+"
        phase_C.Je = - J / 2

        m19_29g = Material("M19_29G")
        m19_29g.conductivity = 1.9e6
        m19_29g.thickness = 0.34
        m19_29g.fill_factor = 0.94
        m19_29g.b = [0.000000, 0.050000, 0.100000, 0.150000, 0.360000, 0.540000, 0.650000, 0.990000, 1.200000, 1.280000,
                     1.330000, 1.360000, 1.440000, 1.520000, 1.580000, 1.630000, 1.670000, 1.800000, 1.900000, 2.000000,
                     2.100000, 2.300000, 2.500000, 2.563994, 3.779890]

        m19_29g.h = [0.000000, 22.280000, 25.460000, 31.830000, 47.740000, 63.660000, 79.570000, 159.150000, 318.300000,
                     477.460000, 636.610000, 795.770000, 1591.500000, 3183.000000, 4774.600000, 6366.100000,
                     7957.700000, 15915.000000, 31830.000000, 111407.000000, 190984.000000, 350135.000000,
                     509252.000000, 560177.200000, 1527756.000000]

        m19_29gsf094 = Material("M19_29GSF094")
        m19_29gsf094.conductivity = 1.9e6
        m19_29gsf094.thickness = 0.34
        m19_29gsf094.fill_factor = 0.94
        m19_29gsf094.b = [0.000000, 0.047002, 0.094002, 0.141002, 0.338404, 0.507605, 0.611006, 0.930612, 1.128024,
                          1.203236, 1.250248, 1.278460, 1.353720, 1.429040, 1.485560, 1.532680, 1.570400, 1.693200,
                          1.788400, 1.888400, 1.988400, 2.188400, 2.388397, 2.452391, 3.668287]

        m19_29gsf094.h = [0.000000, 22.280000, 25.460000, 31.830000, 47.740000, 63.660000, 79.570000, 159.150000,
                          318.300000, 477.460000, 636.610000, 795.770000, 1591.500000, 3183.000000, 4774.600000,
                          6366.100000, 7957.700000, 15915.000000, 31830.000000, 111407.000000, 190984.000000,
                          350135.000000, 509252.000000, 560177.200000, 1527756.000000]

        n36z_50_right = Material("N36Z_50_r")
        n36z_50_right.mu_r = 1.03
        n36z_50_right.coercivity = 782000
        n36z_50_right.conductivity = 0.667e6
        n36z_50_right.remanence_angle = 90 + 90 - self.R3 / 2

        n36z_50_left = copy(n36z_50_right)
        n36z_50_left.name = "N36Z_50_l"
        n36z_50_left.remanence_angle = - n36z_50_right.remanence_angle + 180

        self.snapshot.add_material(air)
        self.snapshot.add_material(phase_A)
        self.snapshot.add_material(phase_B)
        self.snapshot.add_material(phase_C)
        self.snapshot.add_material(m19_29g)
        self.snapshot.add_material(m19_29gsf094)
        self.snapshot.add_material(n36z_50_right)
        self.snapshot.add_material(n36z_50_left)
        self.snapshot.add_material(core)
        self.snapshot.add_material(magnet)

    def assign_materials(self):
        while self.label_queue:
            self.snapshot.assign_material(*self.label_queue.pop())

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)
        pb1 = AntiPeriodicBoundaryCondition("PB1", field_type="magnetic")
        pb2 = AntiPeriodicBoundaryCondition("PB2", field_type="magnetic")
        pb3 = AntiPeriodicBoundaryCondition("PB3", field_type="magnetic")
        pb4 = AntiPeriodicBoundaryCondition("PB4", field_type="magnetic")
        apb = AntiPeriodicAirGap("APairgap", field_type="magnetic", outer_angle=self.rotorangle)

        self.snapshot.add_boundary_condition(a0)
        self.snapshot.add_boundary_condition(pb1)
        self.snapshot.add_boundary_condition(pb2)
        self.snapshot.add_boundary_condition(pb3)
        self.snapshot.add_boundary_condition(pb4)
        self.snapshot.add_boundary_condition(apb)

    def assign_boundary_conditions(self):
        while self.boundary_queue:
            self.snapshot.assign_boundary_condition(*self.boundary_queue.pop())

        while self.boundary_arc_queue:
            self.snapshot.assign_arc_boundary_condition(*self.boundary_arc_queue.pop())

    def _add_slice(self, r_outer, r_inner):
        origin = Node(0.0, 0.0)
        u_left = Node(cos(pi / 2 + 1 / 8 * pi), sin(pi / 2 + 1 / 8 * pi))
        u_right = Node(cos(pi / 2 - 1 / 8 * pi), sin(pi / 2 - 1 / 8 * pi))

        D = u_left * r_outer
        C = u_right * r_outer
        A = u_left * r_inner
        B = u_right * r_inner
        self.add_line(*A, *D)
        self.add_line(*B, *C)
        self.geom.add_arc(CircleArc(C, origin, D, max_seg_deg=1))
        self.geom.add_arc(CircleArc(B, origin, A, max_seg_deg=1))
        return u_left, u_right

    def _add_slices(self):
        # Stator slice
        r_outer = self.S1
        r_inner = self.S2
        ul, ur = self._add_slice(r_outer, r_inner)
        self.label_queue.append((0, self.S1 - 10, "M19_29GSF094"))
        self.boundary_queue.append((*ul * (r_outer + r_inner) / 2, "PB1"))
        self.boundary_queue.append((*ur * (r_outer + r_inner) / 2, "PB1"))

        # stator self.airgap/2 slice
        r_outer = self.S2
        r_inner = self.S2 - self.airgap / 2
        ul, ur = self._add_slice(r_outer, r_inner)
        self.label_queue.append((0, self.S2 - self.airgap / 4, "air"))
        self.boundary_queue.append((*ul * (r_outer + r_inner) / 2, "PB2"))
        self.boundary_queue.append((*ur * (r_outer + r_inner) / 2, "PB2"))
        self.boundary_arc_queue.append((0, r_inner, "APairgap"))

        # rotor slice
        r_outer = self.R1
        r_inner = self.R7
        ul, ur = self._add_slice(r_outer, r_inner)
        self.label_queue.append((0, self.S2 - 5, "M19_29GSF094"))
        self.boundary_queue.append((*ul * (r_outer + r_inner) / 2, "PB4"))
        self.boundary_queue.append((*ur * (r_outer + r_inner) / 2, "PB4"))

        # rotor self.airgap/2 slice
        r_outer = self.R1 + self.airgap / 2
        r_inner = self.R1
        ul, ur = self._add_slice(r_outer, r_inner)
        self.label_queue.append((0, self.R1 + self.airgap / 4, "air"))
        self.boundary_queue.append((*ul * (r_outer + r_inner) / 2, "PB3"))
        self.boundary_queue.append((*ur * (r_outer + r_inner) / 2, "PB3"))
        self.boundary_arc_queue.append((0, r_outer, "APairgap"))

        self.boundary_arc_queue.append((0, self.S1, "a0"))
        self.boundary_arc_queue.append((0, self.R7, "a0"))

    def _add_magnets(self):
        # Magnet Right
        magnet_right = ModelPiece('magnet_right')
        magnet_right.load_piece_from_dxf(self.dir_resources / "clamp.dxf")
        magnet_right.put(self.R5, self.R4, bbox_ref='left')
        magnet_right.rotate(ref_point=magnet_right.left, alpha=-20)

        rect_1 = Rectangle(width=self.R5, height=self.R4)
        magnet_right.geom.add_rectangle(rect_1)
        a = tan((90 - self.R3 / 2) * pi / 180) * (self.R4 - 0.5)
        magnet_right.geom.add_line(Line(Node(0.0, 0.5), Node(-a, 0.5)))
        magnet_right.rotate(ref_point=rect_1.d, alpha=90 - self.R3 / 2)
        magnet_right.translate(dx=0, dy=self.R6 - self.R4)

        magnet_left = copy(magnet_right)
        magnet_left.mirror()
        self.geom.merge_geometry(magnet_right.geom)
        self.geom.merge_geometry(magnet_left.geom)

        self.label_queue.append((10, self.R6, 'N36Z_50_r'))
        self.label_queue.append((-10, self.R6, 'N36Z_50_l'))
        self.label_queue.append((0, 65, 'air'))
        self.label_queue.append((-20, 75, 'air'))
        self.label_queue.append((20, 75, 'air'))

    def _add_slits(self):
        slit = ModelPiece('slit')
        slit.load_piece_from_dxf(self.dir_resources / "slit.dxf")

        # aligning on the y axis
        slit.put(0, 0, bbox_ref='upper')
        slit.rotate(alpha=3.749587282763382)
        slit.put(0, 33.5358, bbox_ref='upper')
        slit.geom.add_line(Line(Node(self.S4 / 2, 0), Node(self.S4 / 2, 1.02579)))
        slit.geom.add_line(Line(Node(-self.S4 / 2, 0), Node(-self.S4 / 2, 1.02579)))

        slit.translate(0, -slit.lower.y + self.S2 - 6.5e-3 + 1.02579)
        slit.rotate(alpha=45 / 2 - 3.75)
        label1 = Node((slit.bbox[0] + slit.bbox[2]) / 2, (slit.bbox[1] + slit.bbox[3]) / 2)
        label2 = slit.lower + 1
        self.label_queue.append((*label1, "A+"))
        self.label_queue.append((*label2, "air"))

        labels = ['A+', 'A+', 'B-', 'B-', 'C+', 'C+']
        for i in range(1, 6):
            slit_i = slit.spawn()
            slit_i.rotate(alpha=-i * 7.5)
            self.geom.merge_geometry(slit_i.geom)

            label1 = label1.rotate(math.radians(-7.5))
            label2 = label2.rotate(math.radians(-7.5))
            self.label_queue.append((*label1, labels[i]))
            self.label_queue.append((*label2, "air"))

        self.geom.merge_geometry(slit.geom)

    def build_geometry(self):
        self._add_slices()

        self._add_magnets()

        self._add_slits()

        self.snapshot.add_geometry(self.geom)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
