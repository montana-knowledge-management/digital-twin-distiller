import math
from copy import copy
from math import cos, pi, radians

from digital_twin_distiller import inch2mm
from digital_twin_distiller.boundaries import DirichletBoundaryCondition
from digital_twin_distiller.boundaries import AntiPeriodicBoundaryCondition
from digital_twin_distiller.boundaries import AntiPeriodicAirGap
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.modelpiece import ModelPiece
from digital_twin_distiller.objects import CircleArc, Node, Line
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot

ModelDir.set_base(__file__)

def cart2pol(x: float, y: float):
    rho = math.hypot(x, y)
    phi = math.atan2(y, x)
    return rho, phi


def pol2cart(rho: float, phi: float):
    x = rho * math.cos(math.radians(phi))
    y = rho * math.sin(math.radians(phi))
    return x, y

ORIGIN = Node(0.0, 0.0)

class PriusMotor(BaseModel):
    """docstring for priusmotor"""
    def __init__(self, **kwargs):
        super(PriusMotor, self).__init__(**kwargs)
        self._init_directories()

        # Geometric parameters
        """ source: http://phdengineeringem.blogspot.com/2018/05/toyota-prius-motor-geometry.html"""
        """ source: https://github.com/Eomys/pyleecan/blob/master/Tests/Data/prius_test.dxf"""
        self.Dso = kwargs.get("Dsi", 269)  # Stator outer diameter [mm]
        self.Dsi = kwargs.get("Dso", 161.93)  # Stator inner diameter [mm]

        self.rotorangle = kwargs.get('rotorangle', 0.0)  # The angle of the sliding band [°]

        self.Dri = kwargs.get("Dri", 111.0)  # Rotor inner diameter [mm]
        self.Dro = kwargs.get("Dro", 160.47)  # Rotor outer diameter [mm]
        self.slheight = kwargs.get("slheight", 7.7)  # Slot height from rotor inner diameter apex point [mm]
        self.mangle = kwargs.get("mangle", 145)  # Magnet angle [°]
        self.mheight = kwargs.get("mheight", 6.5)  # Magnet height [mm]
        self.mwidth = kwargs.get("mwidth", 18.9)  # Magnet width [mm]
        self.aslheight = kwargs.get("aslheight", 1.5)  # Represented on geometry plan.
        self.earheight = kwargs.get("earheight", 1.5)  # Represented on geometry plan.
        self.earlenght1x = kwargs.get("earlenght1x", 2.1)  # Represented on geometry plan.
        self.earlenght2x = kwargs.get("earlenght2x", 1.90)  # Represented on geometry plan.
        self.earlenght2y = kwargs.get("earlenght2y", 2.35)  # Represented on geometry plan.
        self.earlenght3y = kwargs.get("earlenght3y", 1.5)  # Represented on geometry plan.
        self.earlenght4 = kwargs.get("earlenght4", 2.2)  # Represented on geometry plan.

        self.R3 = kwargs.get("R3", 145)  # Magnet material point.
        self.R6 = kwargs.get("R6", 70)  # Magnet material point.

        self.airgap = (self.Dsi - self.Dro) / 2  # Airgap [mm]

        # Excitation setup
        I0 = kwargs.get("I0", 0.0)  # Stator current of one phase [A]
        alpha = kwargs.get("alpha", 0.0)  # Offset of the current [°]
        ina = kwargs.get("ina", 90.0)  # Initial alpha to get T=0Nm position [°]

        coil_area = 0.000142795  # area of the slot [m^2]
        Nturns = 9  # turns of the coil in one slot [u.]
        J0 = Nturns * I0 / coil_area
        self.JU = J0 * cos(radians(ina + alpha))
        self.JV = J0 * cos(radians(ina + alpha + 120))
        self.JW = J0 * cos(radians(ina + alpha + 240))

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False
        femm_metadata.depth = inch2mm(3.3)

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def define_materials(self):

        # define default materials
        air = Material("air")
        air.meshsize = 1.0

        wire = Material("19 AWG")
        wire.lamination_type = "magnetwire"
        wire.diameter = 0.912
        wire.conductivity = 58e6
        wire.meshsize = 1.0

        steel = Material("M19_29GSF094")
        steel.conductivity = 1.9e6
        steel.thickness = 0.34
        steel.fill_factor = 0.94
        steel.b = [0.000000, 0.047002, 0.094002, 0.141002, 0.338404, 0.507605,
                0.611006, 0.930612, 1.128024, 1.203236, 1.250248, 1.278460,
                1.353720, 1.429040, 1.485560, 1.532680, 1.570400, 1.693200,
                1.788400, 1.888400, 1.988400, 2.188400, 2.388397, 2.452391,
                3.668287]

        steel.h = [0.000000, 22.28000, 25.46000, 31.83000, 47.74000, 63.66000,
                79.57000, 159.1500, 318.3000, 477.4600, 636.6100, 795.7700,
                1591.500, 3183.000, 4774.600, 6366.100, 7957.700, 15915.00,
                31830.00, 111407.0, 190984.0, 350135.0, 509252.0, 560177.2,
                1527756.0]

        magnet = Material("N36Z_50")
        magnet.meshsize = 1.0
        magnet.mu_r = 1.03
        magnet.coercivity = 782000
        magnet.conductivity = 0.667e6

        ### create concrete materials
        # Airgap material
        airgap = copy(air)
        airgap.name = 'airgap'
        airgap.meshsize = 0.1

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
        steel_stator.meshsize = 1

        # Rotor steel
        steel_rotor = copy(steel)
        steel_rotor.name = 'steel_rotor'
        steel_rotor.meshsize = 0.3

        # Magnet right
        magnet_right = copy(magnet)
        magnet_right.name = 'magnet_right'
        magnet_right.remanence_angle = -90 + 90 - self.R3 / 2

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

    def add_postprocessing(self):
        entities = [
                (0, 60),
                (0, 68),
                (6, 67),
                (-6, 67),
                (21, 73),
                (-21, 73)
                ]
        self.snapshot.add_postprocessing("integration", entities, "Torque")

    def build_stator_shell(self):

        stator_shell = ModelPiece('stator_shell')

        dsil = Node(*pol2cart(self.Dso / 2, 112.5))
        dsir = Node(*pol2cart(self.Dso / 2, 67.5))

        stator_shell.geom.add_arc(CircleArc(dsir, ORIGIN, dsil, max_seg_deg=10))

        dsol = Node(*pol2cart(self.Dsi / 2, 112.5))
        dsor = Node(*pol2cart(self.Dsi / 2, 67.5))

        stator_shell.geom.add_arc(CircleArc(dsor, ORIGIN, dsol, max_seg_deg=1))

        stator_shell.geom.add_line(Line(dsil, dsol))
        stator_shell.geom.add_line(Line(dsir, dsor))

        agsl = Node(*pol2cart(self.Dsi / 2 - self.airgap / 3, 112.5))
        agsr = Node(*pol2cart(self.Dsi / 2 - self.airgap / 3, 67.5))

        stator_shell.geom.add_arc(CircleArc(agsr, ORIGIN, agsl, max_seg_deg=1))

        stator_shell.geom.add_line(Line(dsol, agsl))
        stator_shell.geom.add_line(Line(dsor, agsr))

        self.geom.merge_geometry(stator_shell.geom)

    def build_slot(self):

        slot= ModelPiece('slot')
        slot.load_piece_from_dxf(ModelDir.RESOURCES / "prius_slot_pyleecan.dxf")
        self.geom.merge_geometry(slot.geom)

    def build_rotor_shell(self):

        rotor_shell = ModelPiece('rotor_shell')

        dril = Node(*pol2cart(self.Dri / 2, 112.5))
        drir = Node(*pol2cart(self.Dri / 2, 67.5))

        rotor_shell.geom.add_arc(CircleArc(drir, ORIGIN, dril, max_seg_deg=10))

        drol = Node(*pol2cart(self.Dro / 2, 112.5))
        dror = Node(*pol2cart(self.Dro / 2, 67.5))

        rotor_shell.geom.add_arc(CircleArc(dror, ORIGIN, drol, max_seg_deg=1))

        rotor_shell.geom.add_line(Line(dril, drol))
        rotor_shell.geom.add_line(Line(drir, dror))

        agrl = Node(*pol2cart(self.Dro / 2 + self.airgap / 3, 112.5))
        agrr = Node(*pol2cart(self.Dro / 2 + self.airgap / 3, 67.5))

        rotor_shell.geom.add_arc(CircleArc(agrr, ORIGIN, agrl, max_seg_deg=1))

        rotor_shell.geom.add_line(Line(drol, agrl))
        rotor_shell.geom.add_line(Line(dror, agrr))

        #rotor_shell.geom.add_node(CircleArc(drir, ORIGIN, dril).apex_pt)

        self.geom.merge_geometry(rotor_shell.geom)

    def build_rotor_slot(self):

        rotor_slot = ModelPiece('rotor_slot')

        temp1 = math.cos(radians((180 - self.mangle) / 2)) * self.mheight

        sorigin = Node(0.0, self.Dri / 2 + self.slheight + temp1)

        pmb1 = Node(0.0, sorigin.y - self.mheight)
        pmb2 = Node(-self.mwidth, pmb1.y)
        pmb3 = Node(-self.mwidth, sorigin.y)

        rotor_slot.geom.add_line(Line(sorigin, pmb1))
        rotor_slot.geom.add_line(Line(pmb2, pmb1))
        rotor_slot.geom.add_line(Line(pmb2, pmb3))
        rotor_slot.geom.add_line(Line(sorigin, pmb3))

        temp2 = math.tan(radians((180 - self.mangle) / 2)) * (self.mheight - self.aslheight)

        apmb1 = Node(0.0, sorigin.y - self.mheight + self.aslheight)
        apmb2 = Node(temp2, apmb1.y)

        rotor_slot.geom.add_line(Line(apmb1, apmb2))

        ear1 = Node(pmb2.x, pmb2.y + self.earheight)
        ear2 = Node(ear1.x - self.earlenght1x, ear1.y)
        ear3 = Node(ear2.x - self.earlenght2x, ear2.y + self.earlenght2y)
        ear4 = Node(ear3.x, ear3.y + self.earlenght3y)
        remy = self.mheight - self.earheight - self.earlenght2y - self.earlenght3y
        remx = math.sqrt((self.earlenght4**2) - (remy**2))
        ear5 = Node(ear4.x + remx, ear4.y + remy)

        rotor_slot.geom.add_line(Line(ear1, ear2))
        rotor_slot.geom.add_line(Line(ear3, ear2))
        rotor_slot.geom.add_line(Line(ear3, ear4))
        rotor_slot.geom.add_line(Line(ear5, ear4))
        rotor_slot.geom.add_line(Line(ear5, pmb3))

        rotor_slot.rotate(ref_point=(sorigin.x, sorigin.y), alpha=-(180 - self.mangle) / 2)

        s = copy(rotor_slot)
        s.mirror(sorigin, ORIGIN)

        #rotor_slot.geom.add_node(sorigin)
        #rotor_slot.geom.add_node(Node(0, sorigin.y - temp1))

        self.geom.merge_geometry(s.geom)
        self.geom.merge_geometry(rotor_slot.geom)

    def build_material(self):
        self.assign_material(10, self.R6, "magnet_right")
        self.assign_material(-10, self.R6, "magnet_left")
        self.assign_material(0, 69, "air")
        self.assign_material(-20, 75, "air")
        self.assign_material(20, 75, "air")

        temp1 = Node(*pol2cart(81.5, 108.75))
        self.assign_material(temp1.x, temp1.y, "air")
        temp2 = Node(*pol2cart(81.5, 101.25))
        self.assign_material(temp2.x, temp2.y, "air")
        temp3 = Node(*pol2cart(81.5, 93.75))
        self.assign_material(temp3.x, temp3.y, "air")
        temp4 = Node(*pol2cart(81.5, 86.25))
        self.assign_material(temp4.x, temp4.y, "air")
        temp5 = Node(*pol2cart(81.5, 78.75))
        self.assign_material(temp5.x, temp5.y, "air")
        temp6 = Node(*pol2cart(81.5, 71.25))
        self.assign_material(temp6.x, temp6.y, "air")

        self.assign_material(0, 79, "steel_rotor")
        self.assign_material(0, 80.35, "air")
        self.assign_material(0, 80.85, "air")
        self.assign_material(0, 120, "steel_stator")

        self.snapshot.add_geometry(self.geom)

    def build_coil(self):

        labels = ["V-", "V-", "U+", "U+", "W-", "W-"]
        label = Node.from_polar(100.0, 71.0)
        for i in range(6):
            self.assign_material(label.x, label.y, labels[i])
            label = label.rotate(pi / 4 / 6)

        self.snapshot.add_geometry(self.geom)

    def build_boundary(self):

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

        self.snapshot.add_geometry(self.geom)

    def build_geometry(self):

        param = 1
        if param == 0:
            rotor = ModelPiece('rotor')
            rotor.load_piece_from_dxf(ModelDir.RESOURCES / "prius_rotor_pyleecan.dxf")
            self.geom.merge_geometry(rotor.geom)

            stator = ModelPiece('stator')
            stator.load_piece_from_dxf(ModelDir.RESOURCES / "prius_stator_pyleecan.dxf")
            self.geom.merge_geometry(stator.geom)

            a = Node.from_polar(80.4494, 67.5)
            b = Node.from_polar(80.4494, 112.5)
            self.geom.add_arc(CircleArc(a,
                                        Node(0,0),
                                        b))
            self.add_line(-30.691, 74.095, -30.7867, 74.3256)
            self.add_line( 30.691, 74.095, 30.7867, 74.3256)

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

            labels = ["V-", "V-","U+", "U+", "W-", "W-"]
            label = Node.from_polar(100.0, 71.0)
            for i in range(6):
                self.assign_material(label.x, label.y, labels[i])
                label = label.rotate(pi/4/6)

            self.assign_boundary(*Node.from_polar(70, 67.5), "PB1")
            self.assign_boundary(*Node.from_polar(70, 112.5), "PB1")

            self.assign_boundary(*Node.from_polar(80.25, 67.5), "PB2")
            self.assign_boundary(*Node.from_polar(80.25, 112.5), "PB2")

            self.assign_boundary(*Node.from_polar(80.8, 67.5), "PB3")
            self.assign_boundary(*Node.from_polar(80.8, 112.5), "PB3")

            self.assign_boundary(*Node.from_polar(110, 67.5), "PB4")
            self.assign_boundary(*Node.from_polar(110, 112.5), "PB4")

            self.assign_boundary_arc(0, 80.4494, "APairgap")
            self.assign_boundary_arc(0,    80.7, "APairgap")

            self.assign_boundary_arc(0, 134.62, "a0")
            self.assign_boundary_arc(0, 55.3199, "a0")

            self.snapshot.add_geometry(self.geom)

        else:

            self.build_stator_shell()
            self.build_slot()
            self.build_rotor_shell()
            self.build_rotor_slot()
            self.build_material()
            self.build_coil()
            self.build_boundary()
            self.snapshot.add_geometry(self.geom)

if __name__ == "__main__":
    m = PriusMotor(exportname="dev")
    print(m(cleanup=True, devmode=True))
