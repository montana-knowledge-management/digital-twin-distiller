from copy import copy

from adze_modeler.model import BaseModel
from adze_modeler.metadata import FemmMetadata
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.material import Material
from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.modelpaths import ModelDir
from adze_modeler.objects import Rectangle
from statistics import fmean

"""
.: https://www.mdpi.com/2076-3417/10/4/1361

.: https://sci-hub.se/10.1109/TPAS.1973.293773

.: https://github.com/tamasorosz/utopya/blob/main/src/optimization_functions_003.py

.: https://sci-hub.se/10.3233/JAE-209504

"""

ModelDir.set_base(__file__)

class PowerTransformer(BaseModel):
    """docstring for PowerTransformer"""
    def __init__(self, **kwargs):
        super(PowerTransformer, self).__init__(**kwargs)
        self._init_directories()

        # Outer rectangle
        self.w1 = kwargs.get('w1', 152)
        self.h1 = kwargs.get('h1', 979+150)
        self.r1 = kwargs.get('r1', 184)
        self.z1 = kwargs.get('z1', 0)

        # LV rectangle
        self.w2 = kwargs.get('w2', 42)
        self.h2 = kwargs.get('h2', 979)
        self.r2 = kwargs.get('r2', 198)
        self.z2 = kwargs.get('z2', 70) # ? 

        # HV rectangle
        self.w3 = kwargs.get('w3', 41)
        self.h3 = kwargs.get('h3', 979)
        self.r3 = kwargs.get('r3', 268)
        self.z3 = kwargs.get('z3', 70)


        # Excitation
        self.jp = kwargs.get('jp', 0.0) * 1_000_000
        self.js = kwargs.get('js', 0.0) * 1_000_000
    
    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "axisymmetric"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = True

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        air = Material('air')
        coil = Material('coil')

        LV = copy(coil)
        LV.name = 'LV'
        LV.Je = self.js

        HV = copy(coil)
        HV.name = 'HV'
        HV.Je = self.jp

        self.snapshot.add_material(air)
        self.snapshot.add_material(LV)
        self.snapshot.add_material(HV)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

    def add_postprocessing(self):
        points = [((self.r1 + self.r2) / 2, self.z1 + self.h1 / 2)]
        self.snapshot.add_postprocessing("integration", points, "Energy")
        
    def build_geometry(self):
        r1 = Rectangle(self.r1, self.z1, width=self.w1, height=self.h1)
        r2 = Rectangle(self.r2, self.z2, width=self.w2, height=self.h2)
        r3 = Rectangle(self.r3, self.z3, width=self.w3, height=self.h3)

        self.geom.add_rectangle(r1)
        self.geom.add_rectangle(r2)
        self.geom.add_rectangle(r3)

        self.snapshot.add_geometry(self.geom)

        self.assign_material(*r2.cp, 'LV')
        self.assign_material(*r3.cp, 'HV')
        self.assign_material((self.r1 + self.r2) / 2, self.z1 + self.h1 / 2, 'air')


        self.assign_boundary(*r1.a.mean(r1.b), 'a0')
        self.assign_boundary(*r1.d.mean(r1.c), 'a0')
        self.assign_boundary(*r1.a.mean(r1.d), 'a0')
        self.assign_boundary(*r1.b.mean(r1.c), 'a0')

    def __call__(self, cleanup=True, devmode=False, timeout=1e5):
        res = super().__call__(cleanup=cleanup, devmode=devmode, timeout=timeout)
        res['ALV'] = self.w2 * self.h2
        res['AHV'] = self.w3 * self.h3
        return res

        

if __name__ == "__main__":
    m = PowerTransformer(exportname="dev", jp=3.0, js=3.02)
    print(m(cleanup=False, devmode=False))
