from copy import copy
import math
from time import perf_counter

from adze_modeler.boundaries import AntiPeriodicAirGap
from adze_modeler.boundaries import AntiPeriodicBoundaryCondition
from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.geometry import Geometry
from adze_modeler.material import Material
from adze_modeler.metadata import FemmMetadata
from adze_modeler.model import BaseModel
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.objects import CircleArc, Line, Node
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.utils import mirror_point

__all__ = ['BLDCMotor']

def cart2pol(x: float, y: float):
    rho = math.hypot(x, y)
    phi = math.atan2(y, x)
    return rho, phi

def pol2cart(rho: float, phi: float):
    x = rho * math.cos(math.radians(phi))
    y = rho * math.sin(math.radians(phi))
    return x, y


ORIGIN = Node(0.0, 0.0)
Y = Node(0.0, 1.0)


class BLDCMotor(BaseModel):
    """
    https://www.femm.info/wiki/RotorMotion

    References:
        .: Ferreira da Luz, M. V., Dular, P., Sadowski, N., Geuzaine, C., & Bastos,
        J. P. A. (2002). Analysis of a permanent magnet generator with dual
        formulations using periodicity conditions and moving band. IEEE
        Transactions on Magnetics, 38(2), 961–964. doi:10.1109/20.996247

        .: Antunes, O. J., Bastos, J. P. A., & Sadowski, N. (2004). Using High-Order
        Finite Elements in Problems With Movement. IEEE Transactions on Magnetics,
        40(2), 529–532. doi:10.1109/tmag.2004.825317

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._init_directories()

        self.rotorangle = kwargs.get('rotorangle', 0.0)
        self.alpha = -kwargs.get('alpha', 0.0)

        # GEOMETRY
        self.depth = kwargs.get('depth', 50.0)

        ## AIRGAP
        self.airgap = kwargs.get('airgap',  0.7)
        self.void = kwargs.get('void',  0.3)

        ## ROTOR
        self.r1 = kwargs.get('r1',  22.8 / 2)  # Rotor Inner Radius
        self.r2 = kwargs.get('r2',  50.5 / 2)  # Rotor Iron Outer Radius
        self.r3 = kwargs.get('r3',  55.1 / 2)  # Rotor Outer Radius
        self.r4 = kwargs.get('r4',  self.r3 + (self.airgap - self.void) / 2)  # Rotor + airgap slice

        ### Magnet
        self.mw = kwargs.get('mw',  15.8566)  # Magnet Width

        ## STATOR
        self.s1 = kwargs.get('s1',  self.r3 + self.airgap)  # Stator Inner Radius
        self.s2 = kwargs.get('s2',  self.s1 + 21.75)  # Stator Outer Radius

        ## SLOT
        self.w1 = kwargs.get('w1',  1.52829)
        self.w2 = kwargs.get('w2',  3.68306)
        self.w3 = kwargs.get('w3',  6.8952)
        self.w4 = kwargs.get('w4',  3.6182)

        self.h1 = kwargs.get('h1',  0.7)
        self.h2 = kwargs.get('h2',  0.3707)
        self.h3 = kwargs.get('h3',  12.1993)
        self.h4 = kwargs.get('h4',  1.9487)

        # Excitation
        coil_area = kwargs.get('coil_area', 7.66533e-5)  # m2
        Nturns = kwargs.get('Nturns', 46)
        I0 = kwargs.get('I0', 0.0)
        J0 = Nturns * I0 / coil_area
        self.JU = J0 * math.cos(math.radians(self.alpha))
        self.JV = J0 * math.cos(math.radians(self.alpha + 120))
        self.JW = J0 * math.cos(math.radians(self.alpha + 240))

        # Mesh sizes
        self.msh_smartmesh = kwargs.get('smartmesh', False)
        self.msh_size_stator_steel = kwargs.get('msh_size_stator_steel', 1.2)
        self.msh_size_rotor_steel = kwargs.get('msh_size_rotor_steel', 0.18)
        self.msh_size_coils = kwargs.get('msh_size_coils', 1.0)
        self.msh_size_air = kwargs.get('msh_size_air', 1.0)
        self.msh_size_airgap = kwargs.get('msh_size_airgap', 0.18)
        self.msh_size_magnets = kwargs.get('msh_size_magnets', 0.18)

        # Materials
        self.m_mur = kwargs.get('magnet_mur', 1.11)
        self.m_Hc = kwargs.get('magnet_Hc', 724000)
        self.m_conductivity = kwargs.get('magnet_conductivity', 1.176e6)
        self.m_angle= kwargs.get('magnet_angle', 90.0)

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = self.msh_smartmesh
        femm_metadata.depth = self.depth

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def add_postprocessing(self):
        points = [
            pol2cart((self.r1 + self.r2) / 2, 90),
            pol2cart((self.r2 + self.r3) / 2, 90),
        ]
        self.snapshot.add_postprocessing("integration", points, "Torque")

    def define_materials(self):
        m19 = Material('M-19 Steel')

        coil = Material('coil')
        coil.meshsize = self.msh_size_coils

        air = Material('air')
        air.meshsize = self.msh_size_air

        smco = Material('SmCo 24 MGOe')
        steel1018 = Material('1018 Steel')

        # Used materials
        stator_steel = copy(m19)
        stator_steel.name = 'stator_steel'
        stator_steel.meshsize = self.msh_size_stator_steel
        stator_steel.thickness = 0.635
        stator_steel.fill_factor = 0.98
        stator_steel.conductivity = 1.9e6
        stator_steel.b = [0.000000, 0.050000, 0.100000, 0.150000, 0.200000,
                          0.250000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000,
                          0.550000, 0.600000, 0.650000, 0.700000, 0.750000, 0.800000,
                          0.850000, 0.900000, 0.950000, 1.000000, 1.050000, 1.100000,
                          1.150000, 1.200000, 1.250000, 1.300000, 1.350000, 1.400000,
                          1.450000, 1.500000, 1.550000, 1.600000, 1.650000, 1.700000,
                          1.750000, 1.800000, 1.850000, 1.900000, 1.950000, 2.000000,
                          2.050000, 2.100000, 2.150000, 2.200000, 2.250000, 2.300000]
        stator_steel.h = [0.000000, 15.120714, 22.718292, 27.842733, 31.871434,
                          35.365044, 38.600588, 41.736202, 44.873979, 48.087807,
                          51.437236, 54.975221, 58.752993, 62.823644, 67.245285,
                          72.084406, 77.420100, 83.350021, 89.999612, 97.537353,
                          106.201406, 116.348464, 128.547329, 143.765431, 163.754169,
                          191.868158, 234.833507, 306.509769, 435.255202, 674.911968,
                          1108.325569, 1813.085468, 2801.217421, 4053.653117,
                          5591.106890, 7448.318413, 9708.815670, 12486.931615,
                          16041.483644, 21249.420624, 31313.495878, 53589.446877,
                          88477.484601, 124329.410540, 159968.569300, 197751.604272,
                          234024.751347]

        # Coils
        # PHASE U
        phase_U_positive = copy(coil)
        phase_U_positive.name = "U+"
        phase_U_positive.Je = self.JU

        phase_U_negative = copy(coil)
        phase_U_negative.name = "U-"
        phase_U_negative.Je = -self.JU

        # PHASE V
        phase_V_positive = copy(coil)
        phase_V_positive.name = "V+"
        phase_V_positive.Je = self.JV

        phase_V_negative = copy(coil)
        phase_V_negative.name = "V-"
        phase_V_negative.Je = -self.JV

        # PHASE W
        phase_W_positive = copy(coil)
        phase_W_positive.name = "W+"
        phase_W_positive.Je = self.JW

        phase_W_negative = copy(coil)
        phase_W_negative.name = "W-"
        phase_W_negative.Je = -self.JW

        ## airgap
        airgap = copy(air)
        airgap.name = 'airgap'
        airgap.meshsize = self.msh_size_airgap

        # Magnet
        magnet = copy(smco)
        magnet.name = 'magnet'
        magnet.meshsize = self.msh_size_magnets
        magnet.mu_r = self.m_mur
        magnet.coercivity = self.m_Hc
        magnet.conductivity = self.m_conductivity
        magnet.remanence_angle = self.m_angle

        # Rotor steel
        rotor_steel = copy(steel1018)
        rotor_steel.name = "rotor_steel"
        rotor_steel.meshsize = self.msh_size_rotor_steel
        rotor_steel.conductivity = 5.8e6
        rotor_steel.b = [0.000000, 0.250300, 0.925000, 1.250000, 1.390000,
                         1.525000, 1.710000, 1.870000, 1.955000, 2.020000, 2.110000,
                         2.225000, 2.430000]
        rotor_steel.h = [0.000000, 238.732500, 795.775000, 1591.550000,
                         2387.325000, 3978.875000, 7957.750000, 15915.500000,
                         23873.250000, 39788.750000, 79577.500000, 159155.000000,
                         318310.000000]
        rotor_steel.phi_hmax = 20

        self.snapshot.add_material(stator_steel)
        self.snapshot.add_material(phase_U_positive)
        self.snapshot.add_material(phase_U_negative)
        self.snapshot.add_material(phase_V_positive)
        self.snapshot.add_material(phase_V_negative)
        self.snapshot.add_material(phase_W_positive)
        self.snapshot.add_material(phase_W_negative)
        self.snapshot.add_material(air)
        self.snapshot.add_material(airgap)
        self.snapshot.add_material(magnet)
        self.snapshot.add_material(rotor_steel)

    def define_boundary_conditions(self):
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

    def build_rotor(self):
        g = Geometry()
        alpha1 = math.degrees(math.asin(self.mw * 0.5 / self.r2))
        alpha2 = math.degrees(math.asin(self.mw * 0.5 / self.r3))

        p1l = Node(*pol2cart(self.r1, 90 + 45 / 2))
        p2l = Node(*pol2cart(self.r2, 90 + 45 / 2))
        p3l = Node(*pol2cart(self.r2, 90 + alpha1))
        p4l = Node(*pol2cart(self.r3, 90 + alpha2))
        p5l = Node(*pol2cart(self.r4, 90 + 45 / 2))

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

        arc1 = CircleArc(p1r, ORIGIN, p1l, max_seg_deg=1)
        arc2 = CircleArc(p5r, ORIGIN, p5l, max_seg_deg=1)
        g.add_arc(arc1)
        g.add_arc(CircleArc(p3l, ORIGIN, p2l, max_seg_deg=1))
        g.add_arc(CircleArc(p2r, ORIGIN, p3r, max_seg_deg=1))
        g.add_arc(CircleArc(p4r, ORIGIN, p4l, max_seg_deg=1))
        g.add_arc(arc2)

        self.geom.merge_geometry(g)

        self.assign_material(0, (self.r1 + self.r2) / 2, 'rotor_steel')
        self.assign_material(0, (self.r2 + self.r3) / 2, 'magnet')
        self.assign_material(*((p2l + p4l) / 2), 'airgap')

        self.assign_boundary_arc(*arc1.apex_pt, 'a0')
        self.assign_boundary_arc(*arc2.apex_pt, "APairgap")
        self.assign_boundary(*((p2l + p5l) / 2), 'PB3')
        self.assign_boundary(*((p2r + p5r) / 2), 'PB3')
        self.assign_boundary(*((p1l + p2l) / 2), 'PB4')
        self.assign_boundary(*((p1r + p2r) / 2), 'PB4')

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

        arc1 = CircleArc(q1r, ORIGIN, q1l, max_seg_deg=1)
        arc2 = CircleArc(q2r, ORIGIN, q2l, max_seg_deg=1)
        arc3 = CircleArc(q3r, ORIGIN, q3l, max_seg_deg=1)
        g.add_arc(arc1)
        g.add_arc(arc2)
        g.add_arc(arc3)

        self.geom.merge_geometry(g)

        self.assign_material(0, self.s1 - (self.airgap - self.void) / 4, 'airgap')
        self.assign_material(0, self.s2 - 0.1, "stator_steel")

        self.assign_boundary_arc(*arc1.apex_pt, "APairgap")
        self.assign_boundary_arc(*arc3.apex_pt, 'a0')
        self.assign_boundary(*((q2l + q3l) / 2), 'PB1')
        self.assign_boundary(*((q2r + q3r) / 2), 'PB1')
        self.assign_boundary(*((q1l + q2l) / 2), 'PB2')
        self.assign_boundary(*((q1r + q2r) / 2), 'PB2')

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

        label1 = (c3l + c4r) / 2
        label2 = (c1l + c2r) / 2

        s.rotate(alpha=-unitangle)
        label1 = label1.rotate(math.radians(-unitangle))
        label2 = label2.rotate(math.radians(-unitangle))
        winding = ['U+', 'V-', 'W+', 'U-', 'V+', 'W-', 'U+', 'V-', 'W+', 'U-',
                   'V+', 'W-', 'U+', 'V-', 'W+', 'U-', 'V+', 'W-', 'U+', 'V-',
                   'W+', 'U-', 'V+', 'W-']
        for i in range(3):
            si = copy(s)
            si.rotate(alpha=i * unitangle)
            self.geom.merge_geometry(si.geom)
            self.assign_material(*label1.rotate(math.radians(i * unitangle)), winding[i])
            self.assign_material(*label2.rotate(math.radians(i * unitangle)), 'air')


    def build_geometry(self):
        self.build_rotor()
        self.build_stator()
        self.build_slots()
        self.snapshot.add_geometry(self.geom)

def execute_model(model: BLDCMotor):
    t0 = perf_counter()
    res = model(timeout=2000, cleanup=False)
    t1 = perf_counter()
    try:
        torque = res["Torque"] * 8
    except Exception as e:
        return None
    # print(f"\t{abs(model.rotorangle):.2f} ° - {abs(model.alpha):.2f} °\t {torque:.3f} Nm \t {t1-t0:.2f} s")
    return torque

if __name__ == "__main__":
    m = BLDCMotor(rotorangle=15 / 4, exportname="dev")
    print(execute_model(m))
