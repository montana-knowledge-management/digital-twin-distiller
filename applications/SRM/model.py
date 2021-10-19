from copy import copy
import math
from time import perf_counter

from adze_modeler.model import BaseModel
from adze_modeler.metadata import FemmMetadata
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.material import Material
from adze_modeler.boundaries import DirichletBoundaryCondition, AntiPeriodicBoundaryCondition, AntiPeriodicAirGap
from adze_modeler.modelpaths import ModelDir
from adze_modeler.geometry import Geometry
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.objects import CircleArc, Line, Node

ModelDir.set_base(__file__)

def cart2pol(x: float, y: float):
    rho = math.hypot(x, y)
    phi = math.atan2(y, x)
    return rho, phi

def pol2cart(rho: float, phi: float):
    x = rho * math.cos(math.radians(phi))
    y = rho * math.sin(math.radians(phi))
    return x, y

def quadroots(a, b, c):
    dis = b * b - 4 * a * c
    sqrt_val = math.sqrt(abs(dis))

    if dis > 0:
        if ((-b - sqrt_val) / (2 * a)) > 0:
            return (-b - sqrt_val) / (2 * a)
        elif ((-b + sqrt_val) / (2 * a)) > 0:
            return (-b + sqrt_val) / (2 * a)
        else:
            return None

    elif dis == 0:
        return(-b / (2 * a))

    else:
        print("Complex Roots")
        return None

def triangle_cos_a(a, b, c):
    angle = math.acos((a ** 2 - b ** 2 - c ** 2) / (-2 * b * c))
    return angle

def triangle_cos_b(a, b, c):
    angle = math.acos((b ** 2 - a ** 2 - c ** 2) / (-2 * a * c))
    return angle

def triangle_cos_c(a, b, c):
    angle = math.acos((c ** 2 - b ** 2 - a ** 2) / (-2 * b * a))
    return angle

class SRM(BaseModel):
    """
    https://livettu-my.sharepoint.com/personal/ekandr_ttu_ee/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fekandr%5Fttu%5Fee%2FDocuments%2FTalTech%2DUWB%20collaboration&originalPath=aHR0cHM6Ly9saXZldHR1LW15LnNoYXJlcG9pbnQuY29tLzpmOi9nL3BlcnNvbmFsL2VrYW5kcl90dHVfZWUvRXZncGlSU1BYdlZLcnlIWk5BYTV3bXNCd0d5Tl96dllvWm94RGpnS3BBcHpYZz9ydGltZT10N2NhRFkyTjJVZw

    TalTech-UWB Collaboration
    """

    def __init__(self, **kwargs):
        exportname = kwargs.get('exportname', None)

        super(SRM, self).__init__()
        self._init_directories()

        ## SIMULATION
        self.rotorangle = kwargs.get('rotorangle',0.0)
        self.alpha = -kwargs.get('alpha', 0.0)
        self.origin = Node(0, 0)

        ## GEOMETRY
        self.depth = kwargs.get('depth', 0.0)
        self.temp1 = kwargs.get('temp1', None)

        ## BOUNDARIES
        self.ag1l = kwargs.get('ag1l', None)
        self.ag1r = kwargs.get('ag1r', None)
        self.ag2l = kwargs.get('ag2l', None)
        self.ag2r = kwargs.get('ag2r', None)

        ## STATOR
        self.D1 = kwargs.get('D1', 103.0)  # Stator outer diameter [mm]
        self.D2 = kwargs.get('D2', 60.0)  # Stator inner diameter [mm]
        self.beta_s = kwargs.get('beta_s', 30.0)  # Stator pole angle [°]
        self.alpha_s = kwargs.get('alpha_s', 90.0)  # Stator tooth angle [°]
        #self.R_bs = kwargs.get('R_bs', 0.0)  # Stator bifurcation radius [mm]
        self.T2 = kwargs.get('T2', 16.5)  # Sloth depth [mm]
        self.gamma_s = kwargs.get('gamma_s', 15.0)  # Coil angle [°]
        self.W1 = kwargs.get('W1', 8.5)  # Coil bottom width [mm]
        self.Z = kwargs.get('Z', 6)  # Number of stator teeth [u.]

        self.p = kwargs.get('p', 2.0)  # Pole pair number

        ## ROTOR
        self.D3 = kwargs.get('D3', 59.5)  # Rotor bore diameter [mm]
        self.D4 = kwargs.get('D4', 10.0)  # Rotor inner diameter [mm]
        self.T1 = kwargs.get('T1', 6.4)  # Core thickness [mm]
        self.beta_r = kwargs.get('beta_r', 29.0)  # Rotor pole angle [°]
        self.alpha_r = kwargs.get('alpha_r', 0.0)  # Rotor tooth angle [°]
        #self.R_br = kwargs.get('R_br', 0.0)  # Rotor bifurcation radius [mm]
        self.gamma_r = kwargs.get('gamma_r', 15.0)  # Rotor pole angle [°]

        ## EXCITATION (Calculated automatically under build_coil)
        self.half_coil_area = kwargs.get('half_coil_area', None)
        self.Nturns = kwargs.get('Nturns', 257.0)
        self.I0 = kwargs.get('I0', 1000.0)

        self.rect1 = kwargs.get('rect1', None)
        self.rect2 = kwargs.get('rect2', None)
        self.rect3 = kwargs.get('rect3', None)
        self.rect4 = kwargs.get('rect4', None)

        ## MESH SIZES
        self.msh_smartmesh = kwargs.get('smartmesh', False)
        self.msh_size_stator_steel = kwargs.get('msh_size_stator_steel', 1.2)
        self.msh_size_rotor_steel = kwargs.get('msh_size_rotor_steel', 0.18)
        self.msh_size_coils = kwargs.get('msh_size_coils', 1.0)
        self.msh_size_air = kwargs.get('msh_size_air', 1.0)
        self.msh_size_airgap = kwargs.get('msh_size_airgap', 0.18)
        self.msh_size_magnets = kwargs.get('msh_size_magnets', 0.18)

        ## MATERIAL


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

    def define_materials(self):  # Coils are defined under build_coil.
        m19 = Material('M-19 Steel')

        self.coil = Material('coil')
        self.coil.meshsize = self.msh_size_coils

        air = Material('air')
        air.meshsize = self.msh_size_air

        ## STATOR STEEL
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

        ## ROTOR STEEL
        rotor_steel = copy(m19)
        rotor_steel.name = 'rotor_steel'
        rotor_steel.meshsize = self.msh_size_rotor_steel
        rotor_steel.thickness = 0.635
        rotor_steel.fill_factor = 0.98
        rotor_steel.conductivity = 1.9e6
        rotor_steel.b = [0.000000, 0.050000, 0.100000, 0.150000, 0.200000,
                          0.250000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000,
                          0.550000, 0.600000, 0.650000, 0.700000, 0.750000, 0.800000,
                          0.850000, 0.900000, 0.950000, 1.000000, 1.050000, 1.100000,
                          1.150000, 1.200000, 1.250000, 1.300000, 1.350000, 1.400000,
                          1.450000, 1.500000, 1.550000, 1.600000, 1.650000, 1.700000,
                          1.750000, 1.800000, 1.850000, 1.900000, 1.950000, 2.000000,
                          2.050000, 2.100000, 2.150000, 2.200000, 2.250000, 2.300000]
        rotor_steel.h = [0.000000, 15.120714, 22.718292, 27.842733, 31.871434,
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

        ## AIRGAP
        airgap = copy(air)
        airgap.name = 'airgap'
        airgap.meshsize = self.msh_size_airgap

        self.snapshot.add_material(stator_steel)
        self.snapshot.add_material(air)
        self.snapshot.add_material(airgap)
        self.snapshot.add_material(rotor_steel)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("A0", field_type="magnetic", magnetic_potential=0.0)
        apb1 = AntiPeriodicBoundaryCondition("APB1", field_type="magnetic")
        apb2 = AntiPeriodicBoundaryCondition("APB2", field_type="magnetic")
        apb3 = AntiPeriodicBoundaryCondition("APB3", field_type="magnetic")
        apb4 = AntiPeriodicBoundaryCondition("APB4", field_type="magnetic")
        apb5 = AntiPeriodicBoundaryCondition("APB5", field_type="magnetic")
        apb6 = AntiPeriodicBoundaryCondition("APB6", field_type="magnetic")
        apb7 = AntiPeriodicBoundaryCondition("APB7", field_type="magnetic")
        slidingband = AntiPeriodicAirGap("slidingband", field_type="magnetic", outer_angle=self.rotorangle)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)
        self.snapshot.add_boundary_condition(apb1)
        self.snapshot.add_boundary_condition(apb2)
        self.snapshot.add_boundary_condition(apb3)
        self.snapshot.add_boundary_condition(apb4)
        self.snapshot.add_boundary_condition(apb5)
        self.snapshot.add_boundary_condition(apb6)
        self.snapshot.add_boundary_condition(apb7)

        self.snapshot.add_boundary_condition(slidingband)

    def add_postprocessing(self):
        half_rotor = (self.D3 / 2 - self.D4 / 2 - self.T1) / 2 + self.D4 / 2 + self.T1

        points = [
            pol2cart((self.D4/2 + self.T1) / 2, 90),
            pol2cart(half_rotor, 90),
            pol2cart(half_rotor, 90 + (360 / (2 * self.p)) / 3),
            pol2cart(half_rotor, 90 - (360 / (2 * self.p)) / 3),
        ]
        self.snapshot.add_postprocessing("integration", points, "Torque")

    def build_rotor(self):
        s = ModelPiece("rotor")

        nr1l = Node(*pol2cart(self.D4/2, 90+(360/(2*self.p))/2))
        nr1r = Node(*pol2cart(self.D4/2, 90-(360/(2*self.p))/2))
        nr2l = Node(*pol2cart(self.D4/2+self.T1, 90+(360/(2*self.p))/2))
        nr2r = Node(*pol2cart(self.D4/2+self.T1, 90-(360/(2*self.p))/2))
        nr3l = Node(*pol2cart(self.D3/2, 90+(360/(2*self.p))/2))
        nr3r = Node(*pol2cart(self.D3/2, 90-(360/(2*self.p))/2))

        s.geom.add_arc(CircleArc(nr1r, self.origin, nr1l, max_seg_deg=1))
        s.geom.add_arc(CircleArc(nr2r, self.origin, nr2l, max_seg_deg=1))
        s.geom.add_arc(CircleArc(nr3r, self.origin, nr3l, max_seg_deg=1))

        s.geom.add_line(Line(nr1l, nr2l))
        s.geom.add_line(Line(nr1r, nr2r))
        s.geom.add_line(Line(nr2l, nr3l))
        s.geom.add_line(Line(nr2r, nr3r))

        self.ag1l = nr3l
        self.ag1r = nr3r

        si = copy(s)
        self.geom.merge_geometry(si.geom)

        self.assign_boundary_arc(*Node(0, self.D4/2), "A0")
        self.assign_boundary(*((nr1l + nr2l) / 2), "APB1")
        self.assign_boundary(*((nr1r + nr2r) / 2), "APB1")
        self.assign_boundary(*((nr2l + nr3l) / 2), "APB2")
        self.assign_boundary(*((nr2r + nr3r) / 2), "APB2")

    def build_rotorpole(self):
        s = ModelPiece("rotorpole")

        nrp2l = Node(*pol2cart(self.D3 / 2, 90 + self.beta_r / 2))
        nrp2r = Node(*pol2cart(self.D3 / 2, 90 - self.beta_r / 2))

        if self.gamma_r < self.beta_r/2:
            a = 1
            b = - 2 * (self.D3 / 2) * math.cos(math.radians(self.beta_r/2 - self.gamma_r))
            c = - ((self.D4 / 2 + self.T1) ** 2 - (self.D3/2) ** 2)
            temp_a = quadroots(a, b, c)

            a = temp_a
            b = self.D4 / 2 + self.T1
            c = self.D3 / 2
            zeta = math.degrees(triangle_cos_a(a, b, c))

            omega = self.beta_r/2 + zeta

        elif self.gamma_r > self.beta_r/2:
            a = 1
            b = - 2 * (self.D3 / 2) * math.cos(math.radians(self.gamma_r - self.beta_r/2))
            c = - ((self.D4 / 2 + self.T1) ** 2 - (self.D3/2) ** 2)
            temp_a = quadroots(a, b, c)

            a = temp_a
            b = self.D4 / 2 + self.T1
            c = self.D3 / 2
            zeta = math.degrees(triangle_cos_a(a, b, c))

            omega = self.beta_r/2 - zeta

        else:
            omega = self.beta_r/2

        nrp1l = Node(*pol2cart(self.D4 / 2 + self.T1, 90 + omega))
        nrp1r = Node(*pol2cart(self.D4 / 2 + self.T1, 90 - omega))

        s.geom.add_line(Line(nrp1r, nrp2r))
        s.geom.add_line(Line(nrp1l, nrp2l))

        si = copy(s)
        self.geom.merge_geometry(si.geom)

    def build_stator(self):
        s = ModelPiece("stator")

        ns1l = Node(*pol2cart(self.D2 / 2, 90 + (360 / (2 * self.p))/2))
        ns1r = Node(*pol2cart(self.D2 / 2, 90 - (360 / (2 * self.p))/2))
        ns2l = Node(*pol2cart(self.D2/2 + self.T2, 90 + (360 / (2 * self.p))/2))
        ns2r = Node(*pol2cart(self.D2/2 + self.T2, 90 - (360 / (2 * self.p))/2))
        ns3l = Node(*pol2cart(self.D1 / 2, 90 + (360 / (2 * self.p))/2))
        ns3r = Node(*pol2cart(self.D1 / 2, 90 - (360 / (2 * self.p))/2))

        s.geom.add_arc(CircleArc(ns1r, self.origin, ns1l, max_seg_deg=1))
        s.geom.add_arc(CircleArc(ns2r, self.origin, ns2l, max_seg_deg=1))
        s.geom.add_arc(CircleArc(ns3r, self.origin, ns3l, max_seg_deg=1))

        s.geom.add_line(Line(ns1l, ns2l))
        s.geom.add_line(Line(ns1r, ns2r))
        s.geom.add_line(Line(ns2l, ns3l))
        s.geom.add_line(Line(ns2r, ns3r))

        self.ag2l = ns1l
        self.ag2r = ns1r

        si = copy(s)
        self.geom.merge_geometry(si.geom)

        self.assign_boundary_arc(*Node(0, self.D1 / 2), "A0")
        self.assign_boundary(*((ns1l + ns2l) / 2), "APB6")
        self.assign_boundary(*((ns1r + ns2r) / 2), "APB6")
        self.assign_boundary(*((ns2l + ns3l) / 2), "APB7")
        self.assign_boundary(*((ns2r + ns3r) / 2), "APB7")


    def build_tooth(self):
        s = ModelPiece("tooth")

        ntl1l = Node(*pol2cart(self.D2 / 2, 90 + self.beta_s/2))
        ntl1r = Node(*pol2cart(self.D2 / 2, 90 - self.beta_s/2))

        if (90-self.alpha_s) < self.beta_s/2:
            a = 1
            b = - 2 * self.D2 / 2 * math.cos(math.radians(270 - self.beta_s / 2 - self.alpha_s))
            c = - ((self.D2 / 2 + self.T2) ** 2 - (self.D2 / 2) ** 2)
            temp_a = quadroots(a, b, c)

            a = temp_a
            b = self.D2 / 2
            c = self.D2 / 2 + self.T2
            zeta = math.degrees(triangle_cos_a(a, b, c))

            omega = self.beta_s / 2 - zeta

        elif (90-self.alpha_s) > self.beta_s/2:
            a = 1
            b = - 2 * self.D2 / 2 * math.cos(math.radians(90 + self.beta_s / 2 + self.alpha_s))
            c = - ((self.D2 / 2 + self.T2) ** 2 - (self.D2 / 2) ** 2)
            temp_a = quadroots(a, b, c)

            a = temp_a
            b = self.D2 / 2
            c = self.D2 / 2 + self.T2
            zeta = math.degrees(triangle_cos_a(a, b, c))

            omega = self.beta_s / 2 + zeta

        else:
            omega = self.beta_s/2

        ntl2l = Node(*pol2cart(self.D2 / 2 + self.T2, 90 + omega))
        ntl2r = Node(*pol2cart(self.D2 / 2 + self.T2, 90 - omega))

        self.temp1 = ntl2r

        s.geom.add_line(Line(ntl1r, ntl2r))
        s.geom.add_line(Line(ntl1l, ntl2l))

        si = copy(s)
        self.geom.merge_geometry(si.geom)

        self.rect1 = ntl1r
        self.rect2 = ntl2r

    def build_coil(self):

        temp = Node(0, self.D2/2 + self.T2)
        x = temp.distance_to(self.temp1)
        theta = math.degrees(math.atan((self.W1+x)/(self.D2/2 + self.T2)))

        ncl2l = Node(*pol2cart(self.D2 / 2 + self.T2, 90 + theta))
        ncl2r = Node(*pol2cart(self.D2 / 2 + self.T2, 90 - theta))

        if self.gamma_s < theta:
            a = 1
            b = - 2 * (self.D2 / 2 + self.T2) * math.cos(math.radians(theta - self.gamma_s))
            c = - ((self.D2 / 2) ** 2 - (self.D2 / 2 + self.T2) ** 2)
            temp_a = quadroots(a, b, c)

            a = temp_a
            b = self.D2 / 2
            c = self.D2 / 2 + self.T2
            zeta = math.degrees(triangle_cos_a(a, b, c))

            omega = theta + zeta

        elif self.gamma_s > theta:
            a = 1
            b = - 2 * (self.D2 / 2 + self.T2) * math.cos(math.radians(self.gamma_s - theta))
            c = - ((self.D2 / 2) ** 2 - (self.D2 / 2 + self.T2) ** 2)
            temp_a = quadroots(a, b, c)

            a = temp_a
            b = self.D2 / 2
            c = self.D2 / 2 + self.T2
            zeta = math.degrees(triangle_cos_a(a, b, c))

            omega = theta - zeta

        else:
            omega = theta

        if omega > self.beta_s/2:
            omega = omega

        else:
            omega = self.beta_s/2

        ncl1l = Node(*pol2cart(self.D2 / 2, 90 + omega))
        ncl1r = Node(*pol2cart(self.D2 / 2, 90 - omega))

        sl = ModelPiece("coil_left")
        sl.geom.add_line(Line(ncl1l, ncl2l))
        sr = ModelPiece("coil_right")
        sr.geom.add_line(Line(ncl1r, ncl2r))

        unitangle = 360/self.Z

        for i in range(2):
            sri = copy(sr)
            sri.rotate(alpha=i * unitangle)
            self.geom.merge_geometry(sri.geom)

        for i in range(2):
            sli = copy(sl)
            sli.rotate(alpha=-i * unitangle)
            self.geom.merge_geometry(sli.geom)

        ## EXCITATION
        self.rect3 = ncl1r
        self.rect4 = ncl2r

        e = self.rect1.distance_to(self.rect4)
        f = self.rect2.distance_to(self.rect3)

        el_dir = self.rect4 - self.rect1
        fl_dir = self.rect2 - self.rect3

        ksi = math.acos((el_dir.x * fl_dir.x + el_dir.y * fl_dir.y)/(e * f))

        self.half_coil_area = (e * f * math.sin(ksi)) / 2 / 1000**2
        J0 = self.Nturns * self.I0 / self.half_coil_area
        self.JU = J0 * math.cos(math.radians(self.alpha))
        self.JV = J0 * math.cos(math.radians(self.alpha + 120))
        self.JW = J0 * math.cos(math.radians(self.alpha + 240))

        # Coils
        # PHASE U
        phase_U_positive = copy(self.coil)
        phase_U_positive.name = "U+"
        phase_U_positive.Je = self.JU

        phase_U_negative = copy(self.coil)
        phase_U_negative.name = "U-"
        phase_U_negative.Je = -self.JU

        # PHASE V
        phase_V_positive = copy(self.coil)
        phase_V_positive.name = "V+"
        phase_V_positive.Je = self.JV

        phase_V_negative = copy(self.coil)
        phase_V_negative.name = "V-"
        phase_V_negative.Je = -self.JV

        # PHASE W
        phase_W_positive = copy(self.coil)
        phase_W_positive.name = "W+"
        phase_W_positive.Je = self.JW

        phase_W_negative = copy(self.coil)
        phase_W_negative.name = "W-"
        phase_W_negative.Je = -self.JW

        self.snapshot.add_material(phase_U_positive)
        self.snapshot.add_material(phase_U_negative)
        self.snapshot.add_material(phase_V_positive)
        self.snapshot.add_material(phase_V_negative)
        self.snapshot.add_material(phase_W_positive)
        self.snapshot.add_material(phase_W_negative)

    def build_boundaries(self):
        s = ModelPiece("boundaries")

        self.nslbl1 = (self.ag2l - self.ag1l) / 3 + self.ag1l
        self.nslbr1 = (self.ag2r - self.ag1r) / 3 + self.ag1r
        self.nslbl2 = (self.ag2l - self.ag1l) / -3 + self.ag2l
        self.nslbr2 = (self.ag2r - self.ag1r) / -3 + self.ag2r

        s.geom.add_line((Line(self.ag1l, self.nslbl1)))
        s.geom.add_line((Line(self.ag1r, self.nslbr1)))
        s.geom.add_line((Line(self.nslbl1, self.nslbl2)))
        s.geom.add_line((Line(self.nslbr1, self.nslbr2)))
        s.geom.add_line((Line(self.nslbl2, self.ag2l)))
        s.geom.add_line((Line(self.nslbr2, self.ag2r)))

        self.slb_rotor = CircleArc(self.nslbr1, self.origin, self.nslbl1)
        self.slb_stator = CircleArc(self.nslbr2, self.origin, self.nslbl2)
        s.geom.add_arc(self.slb_rotor)
        s.geom.add_arc(self.slb_stator)

        si = copy(s)
        self.geom.merge_geometry(si.geom)

        self.assign_boundary(*((self.ag2l - self.ag1l) / 3 / 2 + self.ag1l), "APB3")
        self.assign_boundary(*((self.ag2r - self.ag1r) / 3 / 2 + self.ag1r), "APB3")
        self.assign_boundary(*((self.nslbl1 + self.nslbl2) / 2), "APB4")
        self.assign_boundary(*((self.nslbr1 + self.nslbr2) / 2), "APB4")
        self.assign_boundary(*((self.ag2l - self.ag1l) / -3 / 2 + self.ag2l), "APB5")
        self.assign_boundary(*((self.ag2r - self.ag1r) / -3 / 2 + self.ag2r), "APB5")

        self.assign_boundary_arc(0, self.slb_rotor.apex_pt.y, "slidingband")
        self.assign_boundary_arc(0, self.slb_stator.apex_pt.y, "slidingband")

    def build_material(self):

        self.assign_material(0, (self.slb_rotor.apex_pt.y + self.D3/2) / 2, 'airgap')
        self.assign_material(0, (self.slb_stator.apex_pt.y + self.D2/2) / 2, 'airgap')

        self.assign_material(0, (self.D4 + self.T1)/2, 'rotor_steel')  # rotor

        half_rotor = (self.D3/2 - self.D4/2 - self.T1)/2 + self.D4/2 + self.T1
        self.assign_material(0, half_rotor, 'rotor_steel')  # rotor pole

        self.assign_material(*pol2cart(half_rotor, 90 + (360 / (2 * self.p)) / 3), 'air')  # rotor air left
        self.assign_material(*pol2cart(half_rotor, 90 - (360 / (2 * self.p)) / 3), 'air')  # rotor air right

        self.assign_material(0, (self.D1/2 - self.D2/2 - self.T2)/2 + self.D2/2 + self.T2, 'stator_steel')  # stator rose

        half_stator = (self.D1 / 2 - self.D2 / 2 ) / 2 + self.D2 / 2
        self.assign_material(0, half_stator, 'stator_steel')  # tooth

        self.assign_material(*pol2cart(half_stator, 90 + (360 / (2 * self.p)) / 3), 'air')  # stator air left
        self.assign_material(*pol2cart(half_stator, 90 - (360 / (2 * self.p)) / 3), 'air')  # stator air right

        D1 = (self.rect1.x - self.rect4.x) * (self.rect2.y - self.rect3.y)
        D2 = (self.rect1.y - self.rect4.y) * (self.rect2.x - self.rect3.x)
        D = D1 - D2

        ef_int_x1 = (self.rect1.x * self.rect4.y - self.rect1.y * self.rect4.x) * (self.rect2.x - self.rect3.x)
        ef_int_x2 = (self.rect2.x * self.rect3.y - self.rect2.y * self.rect3.x) * (self.rect1.x - self.rect4.x)
        ef_int_x = (ef_int_x1 - ef_int_x2) / D

        ef_int_y1 = (self.rect1.x * self.rect4.y - self.rect1.y * self.rect4.x) * (self.rect2.y - self.rect3.y)
        ef_int_y2 = (self.rect2.x * self.rect3.y - self.rect2.y * self.rect3.x) * (self.rect1.y - self.rect4.y)
        ef_int_y = (ef_int_y1 - ef_int_y2) / D

        sl = ModelPiece("winding_left")
        label_left = Node(-ef_int_x, ef_int_y)
        sl.geom.add_node(label_left)
        sr = ModelPiece("winding_right")
        label_right = Node(ef_int_x, ef_int_y)
        sr.geom.add_node(label_right)

        unitangle = 360 / self.Z
        winding_left = ['W+', 'V+', 'U+', 'W+', 'V+', 'U+']
        winding_right = ['W-', 'V-', 'U-', 'W-', 'V-', 'U-']

        for i in range(2):
            sli = copy(sl)
            sli.rotate(alpha=-i * unitangle)
            self.assign_material(*label_left.rotate(math.radians(-i * unitangle)), winding_left[-i])

        for i in range(2):
            sri = copy(sr)
            sri.rotate(alpha=i * unitangle)
            self.assign_material(*label_right.rotate(math.radians(i * unitangle)), winding_right[i])

        airgap = (self.D3/2 - self.D2/2)/2 + self.D2/2
        self.assign_material(0, airgap, 'airgap')

    def build_geometry(self):
        self.build_rotor()
        self.build_rotorpole()
        self.build_stator()
        self.build_tooth()
        self.build_coil()
        self.build_boundaries()
        self.build_material()
        self.snapshot.add_geometry(self.geom)

def execute_model(model: SRM):
    t0 = perf_counter()
    res = model(timeout=2000, cleanup=False)
    t1 = perf_counter()
    try:
        torque = res["Torque"] * 6
        print(torque)
    except Exception as e:
        return None
    # print(f"\t{abs(model.rotorangle):.2f} ° - {abs(model.alpha):.2f} °\t {torque:.3f} Nm \t {t1-t0:.2f} s")
    return torque

if __name__ == "__main__":
    m = SRM(exportname="dev")
    print(m(devmode=True))
