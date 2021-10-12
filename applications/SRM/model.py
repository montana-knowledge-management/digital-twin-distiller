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

class SRM(BaseModel):
    """
    https://livettu-my.sharepoint.com/personal/ekandr_ttu_ee/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fekandr%5Fttu%5Fee%2FDocuments%2FTalTech%2DUWB%20collaboration&originalPath=aHR0cHM6Ly9saXZldHR1LW15LnNoYXJlcG9pbnQuY29tLzpmOi9nL3BlcnNvbmFsL2VrYW5kcl90dHVfZWUvRXZncGlSU1BYdlZLcnlIWk5BYTV3bXNCd0d5Tl96dllvWm94RGpnS3BBcHpYZz9ydGltZT10N2NhRFkyTjJVZw

    TalTech-UWB Collaboration
    """

    def __init__(self, **kwargs):
        exportname = kwargs.get('exportname', None)

        super(SRM, self).__init__()
        self._init_directories()

        ## STATOR
        self.D1 = kwargs.get('D1', 0.0)  # Stator outer diameter [mm]
        self.D2 = kwargs.get('D2', 0.0)  # Stator inner diameter [mm]
        self.beta_s = kwargs.get('beta_s', 0.0)  # Stator pole angle [°]
        self.alpha_s = kwargs.get('alpha_s', 0.0)  # Stator tooth angle [°]
        self.R_bs = kwargs.get('R_bs', 0.0)  # Stator bifurcation radius [mm]
        self.T2 = kwargs.get('T2', 0.0)  # Sloth depth [mm]
        self.gamma = kwargs.get('gamma', 0.0)  # Coil angle [°]
        self.W1 = kwargs.get('W1', 0.0)  # Coil bottom width [mm]

        ## ROTOR
        self.D3 = kwargs.get('D3', 0.0)  # Rotor bore diameter [mm]
        self.D4 = kwargs.get('D4', 0.0)  # Rotor inner diameter [mm]
        self.T1 = kwargs.get('T1', 0.0)  # Core thickness [mm]
        self.beta_r = kwargs.get('beta_r', 0.0)  # Rotor pole angle [°]
        self.alpha_r = kwargs.get('alpha_r', 0.0)  # Rotor tooth angle [°]
        self.R_br = kwargs.get('R_br', 0.0)  # Rotor bifurcation radius [mm]






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
        points = [(0, 0), (0, 0)]
        self.snapshot.add_postprocessing("integration", points, "Torque")
        
    def build_geometry(self):
        # ...
        self.snapshot.add_geometry(self.geom)
        

if __name__ == "__main__":
    m = SRM(exportname="dev")
    print(m(cleanup=False))
