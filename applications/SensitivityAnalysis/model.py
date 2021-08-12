from copy import copy
from shutil import rmtree

import numpy as np
from math import pi, sin, cos

from adze_modeler.snapshot import Snapshot
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.model import BaseModel
from adze_modeler.metadata import FemmMetadata
from adze_modeler.platforms.femm import Femm
from adze_modeler.material import Material
from adze_modeler.boundaries import AntiPeriodicBoundaryCondition, DirichletBoundaryCondition
from adze_modeler.objects import Node, CircleArc, Rectangle


class PriusMotor(BaseModel):
    """docstring for PriusMotor"""

    def __init__(self, exportname=None):
        super(PriusMotor, self).__init__(exportname=exportname)

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = True
        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def add_postprocessing(self):
        ...

    def define_materials(self):
        air = Material("air")
        core = Material("core")
        magnet = Material("magnet")

        self.snapshot.add_material(air)
        self.snapshot.add_material(core)
        self.snapshot.add_material(magnet)

    def assign_materials(self):
        ...

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)
        pb1 = AntiPeriodicBoundaryCondition("PB1", field_type="magnetic")
        pb2 = AntiPeriodicBoundaryCondition("PB2", field_type="magnetic")
        pb3 = AntiPeriodicBoundaryCondition("PB3", field_type="magnetic")

        self.snapshot.add_boundary_condition(a0)
        self.snapshot.add_boundary_condition(pb1)
        self.snapshot.add_boundary_condition(pb2)
        self.snapshot.add_boundary_condition(pb3)

    def assign_boundary_conditions(self):
        ...
    
    def add_slice(self, r_outer, r_inner):
        origin = Node(0.0, 0.0)
        u_left = Node(cos(pi / 2 + 3 / 16 * pi), sin(pi / 2 + 3 / 16 * pi))
        u_right = Node(cos(pi / 2 - 3 / 16 * pi), sin(pi / 2 - 3 / 16 * pi))

        D = u_left * r_outer
        C = u_right * r_outer
        A = u_left * r_inner
        B = u_right * r_inner
        self.add_line(*A, *D)
        self.add_line(*B, *C)
        self.geom.add_arc(CircleArc(C, origin, D))
        self.geom.add_arc(CircleArc(B, origin, A))

    def build_geometry(self):
        S1 = 100
        S2 = 60
        R1 = 50
        R3 = 135
        R4 = 6.5
        R5 = 18.9
        R6 = 40
        R7 = 20

        airgap = S2 - R1

        
        # Stator slice
        self.add_slice(S1, S2)

        # rotor slice
        self.add_slice(R1, R7)
        
        # stator airgap/2 slice
        self.add_slice(S2, S2 - airgap / 2)

        # rotor airgap/2 slice
        self.add_slice(R1 + airgap / 2, R1)

        # Magnet Right
        magnet_right = ModelPiece('magnet_right')
        magnet_right.load_piece_from_dxf(self.dir_resources / "clamp.dxf")
        magnet_right.put(R5, R4, bbox_ref='left')
        magnet_right.rotate(ref_point=magnet_right.left, alpha=-20)

        rect_1 = Rectangle(width=R5, height=R4)
        magnet_right.geom.add_rectangle(rect_1)
        magnet_right.rotate(ref_point=rect_1.d, alpha=90-R3/2)
        magnet_right.put(0, R6)

        
        magnet_left = copy(magnet_right)
        magnet_left.mirror()

        
        self.geom.merge_geometry(magnet_right.geom)
        self.geom.merge_geometry(magnet_left.geom)
        self.snapshot.add_geometry(self.geom)


if __name__ == "__main__":
    m = PriusMotor()
    print(m(cleanup=False, devmode=True))
