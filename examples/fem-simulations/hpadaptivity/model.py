from numpy import linspace

from digital_twin_distiller.boundaries import DirichletBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot
from digital_twin_distiller import Node, Line

ModelDir.set_base(__file__)


class ThermalShock(BaseModel):
    """docstring for hpadaptivity"""
    def __init__(self, **kwargs):
        super(ThermalShock, self).__init__(**kwargs)
        self._init_directories()
        self.n_gamma1 = kwargs.get('n_gamma1', 3)
        self.n_gamma2 = kwargs.get('n_gamma2', 3)
        self.n_gamma3 = kwargs.get('n_gamma3', 3)
        self.n_gamma4 = kwargs.get('n_gamma4', 3)

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "heat"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = True
        femm_metadata.depth = 1000

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        air = Material('air')

        self.snapshot.add_material(air)

    def define_boundary_conditions(self):

        for i in range(self.n_gamma1):
            gamma_i = DirichletBoundaryCondition(f"gamma_1_{i+1}", field_type="heat", temperature=0.0)
            self.snapshot.add_boundary_condition(gamma_i)

        for i in range(self.n_gamma2):
            gamma_i = DirichletBoundaryCondition(f"gamma_2_{i+1}", field_type="heat", temperature=0.0)
            self.snapshot.add_boundary_condition(gamma_i)

        for i in range(self.n_gamma3):
            gamma_i = DirichletBoundaryCondition(f"gamma_3_{i+1}", field_type="heat", temperature=0.0)
            self.snapshot.add_boundary_condition(gamma_i)

        for i in range(self.n_gamma4):
            gamma_i = DirichletBoundaryCondition(f"gamma_4_{i+1}", field_type="heat", temperature=0.0)
            self.snapshot.add_boundary_condition(gamma_i)

    def add_postprocessing(self):
        self.snapshot.add_postprocessing("mesh_info", None, None)
        self.snapshot.add_postprocessing("saveimage", None, None)

    def build_geometry(self):
        # ...
        a = Node(0, 0)
        b = Node(1, 0)
        c = Node(1, 1)
        d = Node(0, 1)

        l1 = Line(a, b)
        l2 = Line(b, c)
        l3 = Line(c, d)
        l4 = Line(d, a)
        self.assign_material(0.5, 0.5, "air")
        self.geom.add_line(l1)
        self.geom.add_line(l2)
        self.geom.add_line(l3)
        self.geom.add_line(l4)
        

        self.snapshot.add_geometry(self.geom)


if __name__ == "__main__":
    m = ThermalShock(exportname="dev")
    print(m(cleanup=False, devmode=True))
