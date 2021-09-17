from copy import copy
from pathlib import Path

from adze_modeler.boundaries import DirichletBoundaryCondition, AntiPeriodicBoundaryCondition, AntiPeriodicAirGap
from adze_modeler.material import Material
from adze_modeler.metadata import FemmMetadata
from adze_modeler.model import BaseModel
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.objects import Node, CircleArc, Rectangle, Line
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.utils import inch2mm
import math
from time import perf_counter

__all__= ['DIR_BASE', 'DIR_MEDIA', 'DIR_DATA', 'DIR_RESOURCES', 'DIR_SNAPSHOTS', 'PriusMotor',
          'ParametricPriusMotor', 'execute_model']

DIR_BASE = Path(__file__).parent
DIR_MEDIA = DIR_BASE / "media"
DIR_DATA = DIR_BASE / "data"
DIR_RESOURCES = DIR_BASE / "resources"
DIR_SNAPSHOTS = DIR_BASE / "snapshots"


class PriusMotor(BaseModel):
    def __init__(self, rotorangle=0.0, alpha=0.0, I0=0.0, exportname: str = None):
        """
        Parameters:
                rotorangle: Then angle of the rotor
                alpha: Electrical angle
                I0: The magnitude of the excitation current
                exportname: Specific directory name
        """
        super().__init__(exportname=exportname)
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

        # Mechanical and Electrical angles
        self.rotorangle = -rotorangle
        self.alpha = -alpha

        # Excitation setup
        coil_area = 0.000142795  # m2
        Nturns = 9
        J0 = Nturns * I0 / coil_area
        self.JU = J0 * math.cos(math.radians(alpha))
        self.JV = J0 * math.cos(math.radians(alpha + 120))
        self.JW = J0 * math.cos(math.radians(alpha + 240))

        # Mesh sizes
        self.msh_size_stator_steel = 1.2
        self.msh_size_rotor_steel = 0.5
        self.msh_size_coils = 1.0
        self.msh_size_air = 1.0
        self.msh_size_airgap = 0.3
        self.msh_size_magnets = 1.0

        self.offset=0.25

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

    def add_postprocessing(self):
        entities = [(0, self.R6 - 0.1),
                    (0, self.R1 - 0.1),
                    (-10, 70),
                    (10, 70),
                    (-20, 75),
                    (20, 75)]
        self.snapshot.add_postprocessing("integration", entities, "Torque")

    def define_materials(self):
        # define default materials
        air = Material("air")
        air.meshsize = self.msh_size_air

        wire = Material("20 AWG")
        wire.lamination_type = "magnetwire"
        wire.diameter = 0.812049969500513
        wire.conductivity = 58e6
        wire.meshsize = self.msh_size_coils

        steel = Material("M19_29GSF094")
        steel.conductivity = 1.9e6
        steel.thickness = 0.34
        steel.fill_factor = 0.94
        steel.b = [0.000000, 0.047002, 0.094002, 0.141002, 0.338404, 0.507605, 0.611006, 0.930612, 1.128024, 1.203236,
                   1.250248, 1.278460, 1.353720, 1.429040, 1.485560, 1.532680, 1.570400, 1.693200, 1.788400, 1.888400,
                   1.988400, 2.188400, 2.388397, 2.452391, 3.668287]

        steel.h = [0.000000, 22.280000, 25.460000, 31.830000, 47.740000, 63.660000, 79.570000, 159.150000, 318.300000,
                   477.460000, 636.610000, 795.770000, 1591.500000, 3183.000000, 4774.600000, 6366.100000, 7957.700000,
                   15915.000000, 31830.000000, 111407.000000, 190984.000000, 350135.000000, 509252.000000,
                   560177.200000, 1527756.000000]

        magnet = Material("N36Z_50")
        magnet.meshsize = self.msh_size_magnets
        magnet.mu_r = 1.03
        magnet.coercivity = 782000
        magnet.conductivity = 0.667e6

        ### create concrete materials
        # Airgap material
        airgap = copy(air)
        airgap.name = 'airgap'
        airgap.meshsize = self.msh_size_airgap

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
        steel_stator.name = 'steel_stator'
        steel_stator.meshsize = self.msh_size_stator_steel

        # Rotor steel
        steel_rotor = copy(steel)
        steel_rotor.name = 'steel_rotor'
        steel_rotor.meshsize = self.msh_size_stator_steel

        # Magnet right
        magnet_right = copy(magnet)
        magnet_right.name = 'magnet_right'
        magnet_right.remanence_angle = 90 + 90 - self.R3 / 2

        # Magnet left
        magnet_left = copy(magnet)
        magnet_left.name = 'magnet_left'
        magnet_left.remanence_angle = -magnet_right.remanence_angle + 180

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

    def _add_slice(self, r_outer, r_inner, segment_deg=1):
        origin = Node(0.0, 0.0)
        phi_left = math.pi / 2 + 1 / 8 * math.pi
        phi_right = math.pi / 2 - 1 / 8 * math.pi
        u_left = Node(math.cos(phi_left), math.sin(phi_left))
        u_right = Node(math.cos(phi_right), math.sin(phi_right))

        D = (u_left * r_outer)
        C = (u_right * r_outer)
        A = (u_left * r_inner)
        B = (u_right * r_inner)
        self.add_line(*A, *D)
        self.add_line(*B, *C)
        self.geom.add_arc(CircleArc(C, origin, D, max_seg_deg=segment_deg))
        self.geom.add_arc(CircleArc(B, origin, A, max_seg_deg=segment_deg))
        return u_left, u_right

    def _add_slices(self):
        # Stator slice
        r_outer = self.S1+self.offset
        r_inner = self.S2+self.offset
        ul, ur = self._add_slice(r_outer, r_inner)
        label_boudnary_left = ul * (r_outer + r_inner) / 2
        label_boudnary_right = ur * (r_outer + r_inner) / 2
        self.assign_material(0, self.S1 - 10+self.offset, "steel_stator")
        self.assign_boundary(*label_boudnary_left, "PB1")
        self.assign_boundary(*label_boudnary_right, "PB1")

        # stator  - airgap / 2 slice
        r_outer = self.S2+self.offset
        r_inner = self.S2 - self.airgap / 2+self.offset
        ul, ur = self._add_slice(r_outer, r_inner, segment_deg=1)
        label_boudnary_left = ul * (r_outer + r_inner) / 2
        label_boudnary_right = ur * (r_outer + r_inner) / 2
        self.assign_material(0, self.S2 - self.airgap / 4+self.offset, "airgap")
        self.assign_boundary(*label_boudnary_left, "PB2")
        self.assign_boundary(*label_boudnary_right, "PB2")
        self.assign_boundary_arc(0, r_inner, "APairgap")

        # rotor slice
        r_outer = self.R1
        r_inner = self.R7
        ul, ur = self._add_slice(r_outer, r_inner)
        label_boudnary_left = ul * (r_outer + r_inner) / 2
        label_boudnary_right = ur * (r_outer + r_inner) / 2
        self.assign_material(0, self.S2 - 5, "steel_rotor")
        self.assign_boundary(*label_boudnary_left, "PB4")
        self.assign_boundary(*label_boudnary_right, "PB4")

        # rotor + airgap / 2 slice
        r_outer = self.R1 + self.airgap / 2
        r_inner = self.R1
        ul, ur = self._add_slice(r_outer, r_inner, segment_deg=1)
        label_boudnary_left = ul * (r_outer + r_inner) / 2
        label_boudnary_right = ur * (r_outer + r_inner) / 2
        self.assign_material(0, self.R1 + self.airgap / 4, "airgap")
        self.assign_boundary(*label_boudnary_left, "PB3")
        self.assign_boundary(*label_boudnary_right, "PB3")
        self.assign_boundary_arc(0, r_outer, "APairgap")

        self.assign_boundary_arc(0, self.S1, "a0")
        self.assign_boundary_arc(0, self.R7, "a0")

    def _add_magnets(self):
        # Magnet Right
        magnet_right = ModelPiece("magnet_right")
        magnet_right.load_piece_from_dxf(self.dir_resources / "clamp.dxf")
        magnet_right.put(self.R5, self.R4, bbox_ref="left")
        magnet_right.rotate(ref_point=magnet_right.left, alpha=-20)

        rect_1 = Rectangle(width=self.R5, height=self.R4)
        magnet_right.geom.add_rectangle(rect_1)
        a = math.tan(math.radians(90 - self.R3 / 2)) * (self.R4 - 0.5)
        magnet_right.geom.add_line(Line(Node(0.0, 0.5), Node(-a, 0.5)))
        magnet_right.rotate(ref_point=rect_1.d, alpha=90 - self.R3 / 2)
        magnet_right.translate(dx=0, dy=self.R6 - self.R4)

        magnet_left = copy(magnet_right)
        magnet_left.mirror()
        self.geom.merge_geometry(magnet_right.geom)
        self.geom.merge_geometry(magnet_left.geom)

        self.assign_material(10, self.R6, "magnet_right")
        self.assign_material(-10, self.R6, "magnet_left")
        self.assign_material(0, 65, "air")
        self.assign_material(-20, 75, "air")
        self.assign_material(20, 75, "air")

    def _add_slits(self):
        slit = ModelPiece("slit")
        slit.load_piece_from_dxf(self.dir_resources / "slit.dxf")

        # aligning on the y axis
        slit.put(0, 0, bbox_ref="upper")
        slit.rotate(alpha=3.749587282763382)
        slit.put(0, 33.5358, bbox_ref="upper")
        slit.geom.add_line(Line(Node(self.S4 / 2, 0), Node(self.S4 / 2, 1.02579)))
        slit.geom.add_line(Line(Node(-self.S4 / 2, 0), Node(-self.S4 / 2, 1.02579)))

        slit.translate(0, self.S2 - 5.8e-3+self.offset)
        slit.rotate(alpha=45 / 2 - 3.75)
        label1 = Node((slit.bbox[0] + slit.bbox[2]) / 2, (slit.bbox[1] + slit.bbox[3]) / 2)
        label2 = slit.lower + 1
        self.assign_material(*label1, "U+")
        self.assign_material(*label2, "air")

        labels = ["U+", "V-", "V-", "W+", "W+", "U-"]
        for i in range(1, 6):
            slit_i = slit.spawn()
            slit_i.rotate(alpha=-i * 7.5)
            self.geom.merge_geometry(slit_i.geom)

            label1 = label1.rotate(math.radians(-7.5))
            label2 = label2.rotate(math.radians(-7.5))
            self.assign_material(*label1, labels[i])
            self.assign_material(*label2, "air")

        self.geom.merge_geometry(slit.geom)

    def build_geometry(self):
        self._add_slices()

        self._add_magnets()

        self._add_slits()

        self.snapshot.add_geometry(self.geom)

    def __repr__(self):
        return f"{self.rotorangle:.2f} ° - {self.alpha:.2f}°"


class ParametricPriusMotor(PriusMotor):

    def __init__(self, X: dict, rotorangle=0.0, alpha=0.0, I0=0.0, exportname=None):
        super(ParametricPriusMotor, self).__init__(rotorangle=rotorangle,
                                                   alpha=alpha, I0=I0, exportname=exportname)

        self.S1 += X.get("S1", 0.0)
        self.S2 += X.get("S2", 0.0)
        self.S4 += X.get("S4", 0.0)
        self.R1 += X.get("R1", 0.0)
        self.R3 += X.get("R3", 0.0)
        self.R4 += X.get("R4", 0.0)
        self.R5 += X.get("R5", 0.0)
        self.R6 += X.get("R6", 0.0)
        self.R7 += X.get("R7", 0.0)

        self.airgap += X.get("airgap", 0.0)

def execute_model(model: PriusMotor):
    t0 = perf_counter()
    res = model(timeout=2000, cleanup=True)
    t1 = perf_counter()
    torque = res["Torque"]*8
    print(f"\t {abs(model.rotorangle):.2f} ° - {abs(model.alpha):.2f} °\t {torque:.3f} Nm \t {t1-t0:.2f} s")
    return torque

if __name__ == "__main__":
    from numpy import linspace
    import multiprocessing

    m = PriusMotor(rotorangle=360/48*3/4, exportname="dev")
    execute_model(m)