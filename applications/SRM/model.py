from copy import copy
import math
from time import perf_counter

from adze_modeler.model import BaseModel
from adze_modeler.metadata import FemmMetadata
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.material import Material
from adze_modeler.boundaries import DirichletBoundaryCondition
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
        self.rotorangle = kwargs.get('rotorangle', 0.0)
        self.alpha = -kwargs.get('alpha', 0.0)
        self.origin = Node(0, 0)

        ## GEOMETRY
        self.depth = kwargs.get('depth', 0.0)

        ## STATOR
        self.D1 = kwargs.get('D1', 103.0)  # Stator outer diameter [mm]
        self.D2 = kwargs.get('D2', 60.0)  # Stator inner diameter [mm]
        self.beta_s = kwargs.get('beta_s', 30.0)  # Stator pole angle [°]
        self.alpha_s = kwargs.get('alpha_s', 0.0)  # Stator tooth angle [°]
        self.R_bs = kwargs.get('R_bs', 0.0)  # Stator bifurcation radius [mm]
        self.T2 = kwargs.get('T2', 33.0)  # Sloth depth [mm]
        self.gamma = kwargs.get('gamma', 0.0)  # Coil angle [°]
        self.W1 = kwargs.get('W1', 0.0)  # Coil bottom width [mm]

        self.p = kwargs.get('p', 3.0)  # Pole pair number
        ## ROTOR
        self.D3 = kwargs.get('D3', 0.0)  # Rotor bore diameter [mm]
        self.D4 = kwargs.get('D4', 10.0)  # Rotor inner diameter [mm]
        self.T1 = kwargs.get('T1', 14.0)  # Core thickness [mm]
        self.beta_r = kwargs.get('beta_r', 0.0)  # Rotor pole angle [°]
        self.alpha_r = kwargs.get('alpha_r', 0.0)  # Rotor tooth angle [°]
        self.R_br = kwargs.get('R_br', 0.0)  # Rotor bifurcation radius [mm]

        ## EXCITATION
        half_coil_area = kwargs.get('half_coil_area', 1.0)  # m2
        Nturns = kwargs.get('Nturns', 0.0)
        I0 = kwargs.get('I0', 0.0)
        J0 = Nturns * I0 / half_coil_area
        self.JU = J0 * math.cos(math.radians(self.alpha))
        self.JV = J0 * math.cos(math.radians(self.alpha + 120))
        self.JW = J0 * math.cos(math.radians(self.alpha + 240))

        ## MESH SIZES
        self.msh_smartmesh = kwargs.get('smartmesh', False)
        self.msh_size_stator_steel = kwargs.get('msh_size_stator_steel', 0.0)
        self.msh_size_rotor_steel = kwargs.get('msh_size_rotor_steel', 0.0)
        self.msh_size_coils = kwargs.get('msh_size_coils', 0.0)
        self.msh_size_air = kwargs.get('msh_size_air', 0.0)
        self.msh_size_airgap = kwargs.get('msh_size_airgap', 0.0)
        self.msh_size_magnets = kwargs.get('msh_size_magnets', 0.0)

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

    def define_materials(self):
        air = Material('air')

        self.snapshot.add_material(air)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

    def add_postprocessing(self):
        #points = [(0, 0), (0, 0)]
        #self.snapshot.add_postprocessing("integration", points, "Torque")
        pass

    def build_rotor(self):
        s = ModelPiece("rotor")

        nr1l = Node(*pol2cart(self.D4/2, 90+360/(2*self.p)))
        nr1r = Node(*pol2cart(self.D4/2, 90-360/(2*self.p)))
        nr2l = Node(*pol2cart((self.D4+self.T1)/2, 90+360/(2*self.p)))
        nr2r = Node(*pol2cart((self.D4+self.T1)/2, 90-360/(2*self.p)))

        s.geom.add_arc(CircleArc(nr1r, self.origin, nr1l, max_seg_deg=1))
        s.geom.add_arc(CircleArc(nr2r, self.origin, nr2l, max_seg_deg=1))

        si = copy(s)
        self.geom.merge_geometry(si.geom)

    def build_stator(self):
        s = ModelPiece("stator")

        ns1l = Node(*pol2cart(self.D2 / 2, 90 + 360 / (2 * self.p)))
        ns1r = Node(*pol2cart(self.D2 / 2, 90 - 360 / (2 * self.p)))
        ns2l = Node(*pol2cart((self.D2 + self.T2) / 2, 90 + 360 / (2 * self.p)))
        ns2r = Node(*pol2cart((self.D2 + self.T2) / 2, 90 - 360 / (2 * self.p)))
        ns3l = Node(*pol2cart(self.D1 / 2, 90 + 360 / (2 * self.p)))
        ns3r = Node(*pol2cart(self.D1 / 2, 90 - 360 / (2 * self.p)))

        s.geom.add_arc(CircleArc(ns1r, self.origin, ns1l, max_seg_deg=1))
        s.geom.add_arc(CircleArc(ns2r, self.origin, ns2l, max_seg_deg=1))
        s.geom.add_arc(CircleArc(ns3r, self.origin, ns3l, max_seg_deg=1))

        si = copy(s)
        self.geom.merge_geometry(si.geom)

    def build_slot(self):
        s = ModelPiece("stator")

        nsl1l = Node(*pol2cart(self.D2 / 2, 90 + self.beta_s))
        nsl1r = Node(*pol2cart(self.D2 / 2, 90 - self.beta_s))
        #nsl2l =
        #nsl2r =

        #s.geom.add_line(Line(nsl1r, nsl2r))
        #s.geom.add_line(Line(nsl1l, nsl2l))

        si = copy(s)
        self.geom.merge_geometry(si.geom)

    def build_geometry(self):
        self.build_rotor()
        self.build_stator()
        self.build_slot()
        self.snapshot.add_geometry(self.geom)

if __name__ == "__main__":
    m = SRM(exportname="dev")
    print(m(cleanup=False))
