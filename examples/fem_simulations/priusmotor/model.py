from copy import copy
from math import cos, pi, radians

from digital_twin_distiller import Line, inch2mm
from digital_twin_distiller.boundaries import (
    AntiPeriodicAirGap,
    AntiPeriodicBoundaryCondition,
    DirichletBoundaryCondition,
)
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.modelpiece import ModelPiece
from digital_twin_distiller.objects import CircleArc, Node
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot

ModelDir.set_base(__file__)


class PriusMotor(BaseModel):
    """docstring for priusmotor"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._init_directories()

        # Geometric parameters
        self.S1 = inch2mm(10.600) / 2
        self.S2 = inch2mm(6.375) / 2
        self.S4 = inch2mm(0.076)
        self.R1 = inch2mm(6.315) / 2
        self.R3 = 145
        self.R4 = 6.5
        self.R5 = 18.9
        self.R6 = 70.0
        self.R7 = inch2mm(4.356) / 2
        self.airgap = self.S2 - self.R1

        self.rotorangle = kwargs.get("rotorangle", 0.0)

        # Excitation setup
        I0 = kwargs.get("I0", 0.0)
        alpha = kwargs.get("alpha", 0.0)

        coil_area = 0.000142795  # m2
        Nturns = 9
        J0 = Nturns * I0 / coil_area
        self.JU = J0 * cos(radians(alpha))
        self.JV = J0 * cos(radians(alpha + 120))
        self.JW = J0 * cos(radians(alpha + 240))

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False
        femm_metadata.depth = 83.6

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def define_materials(self):

        # define default materials
        air = Material("air")
        air.meshsize = 1.0

        wire = Material("20 AWG")
        wire.lamination_type = "magnetwire"
        wire.diameter = 0.812049969500513
        wire.conductivity = 58e6
        wire.meshsize = 1.0

        steel = Material("M19_29GSF094")
        steel.conductivity = 1.9e6
        steel.thickness = 0.34
        steel.fill_factor = 0.94
        steel.b = [
            0.000000,
            0.047002,
            0.094002,
            0.141002,
            0.338404,
            0.507605,
            0.611006,
            0.930612,
            1.128024,
            1.203236,
            1.250248,
            1.278460,
            1.353720,
            1.429040,
            1.485560,
            1.532680,
            1.570400,
            1.693200,
            1.788400,
            1.888400,
            1.988400,
            2.188400,
            2.388397,
            2.452391,
            3.668287,
        ]

        steel.h = [
            0.0,
            22.28,
            25.46,
            31.83,
            47.74,
            63.66,
            79.57,
            159.15,
            318.3,
            477.46,
            636.61,
            795.77,
            1591.5,
            3183.0,
            4774.6,
            6366.1,
            7957.7,
            15915.0,
            31830.0,
            111407.000000,
            190984.000000,
            350135.0,
            509252.0,
            560177.2,
            1527756.0,
        ]

        magnet = Material("N36Z_50")
        magnet.meshsize = 1.0
        magnet.mu_r = 1.03
        magnet.coercivity = 782000
        magnet.conductivity = 0.667e6

        ### create concrete materials
        # Airgap material
        airgap = copy(air)
        airgap.name = "airgap"
        airgap.meshsize = 0.3

        # Coils
        # PHASE U
        phase_U_positive = copy(wire)
        phase_U_positive.name = "U+"
        phase_U_positive.Je = self.JU

        phase_U_negative = copy(wire)
        phase_U_negative.name = "U-"
        phase_U_negative.Je = -self.JU

        # PHASE V
        phase_V_positive = copy(wire)
        phase_V_positive.name = "V+"
        phase_V_positive.Je = self.JV

        phase_V_negative = copy(wire)
        phase_V_negative.name = "V-"
        phase_V_negative.Je = -self.JV

        # PHASE W
        phase_W_positive = copy(wire)
        phase_W_positive.name = "W+"
        phase_W_positive.Je = self.JW

        phase_W_negative = copy(wire)
        phase_W_negative.name = "W-"
        phase_W_negative.Je = -self.JW

        # Stator steel
        steel_stator = copy(steel)
        steel_stator.name = "steel_stator"
        steel_stator.meshsize = 1.2

        # Rotor steel
        steel_rotor = copy(steel)
        steel_rotor.name = "steel_rotor"
        steel_rotor.meshsize = 0.5

        # Magnet right
        magnet_right = copy(magnet)
        magnet_right.name = "magnet_right"
        magnet_right.remanence_angle = 107.46

        # Magnet left
        magnet_left = copy(magnet)
        magnet_left.name = "magnet_left"
        magnet_left.remanence_angle = 72.54

        # Adding the used materials to the snapshot
        self.snapshot.add_material(air)
        self.snapshot.add_material(airgap)
        self.snapshot.add_material(phase_U_positive)
        self.snapshot.add_material(phase_U_negative)
        self.snapshot.add_material(phase_V_positive)
        self.snapshot.add_material(phase_V_negative)
        self.snapshot.add_material(phase_W_positive)
        self.snapshot.add_material(phase_W_negative)
        self.snapshot.add_material(steel_stator)
        self.snapshot.add_material(steel_rotor)
        self.snapshot.add_material(magnet_right)
        self.snapshot.add_material(magnet_left)

    def define_boundary_conditions(self):
        # Define boundary conditions
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)
        pb1 = AntiPeriodicBoundaryCondition("PB1", field_type="magnetic")
        pb2 = AntiPeriodicBoundaryCondition("PB2", field_type="magnetic")
        pb3 = AntiPeriodicBoundaryCondition("PB3", field_type="magnetic")
        pb4 = AntiPeriodicBoundaryCondition("PB4", field_type="magnetic")
        apb = AntiPeriodicAirGap("APairgap", field_type="magnetic", outer_angle=self.rotorangle)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)
        self.snapshot.add_boundary_condition(pb1)
        self.snapshot.add_boundary_condition(pb2)
        self.snapshot.add_boundary_condition(pb3)
        self.snapshot.add_boundary_condition(pb4)
        self.snapshot.add_boundary_condition(apb)

    def add_postprocessing(self):
        entities = [(0, 60), (0, 65), (6, 67), (-6, 67), (21, 73), (-21, 73)]
        self.snapshot.add_postprocessing("integration", entities, "Torque")

    def build_geometry(self):
        # ...
        rotor = ModelPiece("rotor")
        rotor.load_piece_from_svg("/home/csanyig/PycharmProjects/bosch-2022-1/resources/prius_motor/rotor.svg")
        # rotor = ModelPiece('rotor')
        # rotor.load_piece_from_dxf(ModelDir.RESOURCES / "prius_test.dxf")
        # rotor.rotate(alpha=67.5)
        # rotor.scale(1000, 1000)
        self.geom.merge_geometry(rotor.geom)

        stator = ModelPiece("stator")
        stator.load_piece_from_svg("/home/csanyig/PycharmProjects/bosch-2022-1/resources/prius_motor/stator.svg")
        # stator.load_piece_from_dxf(ModelDir.RESOURCES / "Prius2004_Stator_SB.dxf")
        self.geom.merge_geometry(stator.geom)

        a = Node.from_polar(80.4494, 67.5)
        b = Node.from_polar(80.4494, 112.5)
        self.geom.add_arc(CircleArc(a, Node(0, 0), b))
        self.add_line(-30.691, 74.095, -30.7867, 74.3256)
        self.add_line(30.691, 74.095, 30.7867, 74.3256)

        a = Node.from_polar(80.8, 67.5)
        b = Node.from_polar(80.8, 112.5)
        self.geom.add_arc(CircleArc(a, Node(0, 0), b))
        self.add_line(-30.9782, 74.788, *b)
        self.add_line(30.9782, 74.788, *a)

        self.snapshot.add_geometry(self.geom)
        for i in range(len(self.geom.circle_arcs)):
            self.geom.circle_arcs[i].max_seg_deg = 1

        self.assign_material(10, self.R6, "magnet_right")
        self.assign_material(-10, self.R6, "magnet_left")
        self.assign_material(0, 65, "air")
        self.assign_material(-20, 75, "air")
        self.assign_material(20, 75, "air")
        self.assign_material(-5.5, 81.3, "air")

        self.assign_material(0, 79, "steel_rotor")
        self.assign_material(0, 80.28, "air")
        self.assign_material(0, 120, "steel_stator")

        labels = ["U+", "V-", "V-", "W+", "W+", "U-"]
        label = Node.from_polar(100.0, 71.0)
        for i in range(6):
            self.assign_material(label.x, label.y, labels[i])
            label = label.rotate(pi / 4 / 6)

        self.assign_boundary(*Node.from_polar(70, 67.5), "PB1")
        self.assign_boundary(*Node.from_polar(70, 112.5), "PB1")

        self.assign_boundary(*Node.from_polar(80.25, 67.5), "PB2")
        self.assign_boundary(*Node.from_polar(80.25, 112.5), "PB2")

        self.assign_boundary(*Node.from_polar(80.8, 67.5), "PB3")
        self.assign_boundary(*Node.from_polar(80.8, 112.5), "PB3")

        self.assign_boundary(*Node.from_polar(110, 67.5), "PB4")
        self.assign_boundary(*Node.from_polar(110, 112.5), "PB4")

        self.assign_boundary_arc(0, 80.4494, "APairgap")
        self.assign_boundary_arc(0, 80.7, "APairgap")

        self.assign_boundary_arc(0, 134.62, "a0")
        self.assign_boundary_arc(0, 55.3199, "a0")


if __name__ == "__main__":
    m = PriusMotor(exportname="dev")
    print(m(cleanup=False, devmode=False, timeout=5000))
