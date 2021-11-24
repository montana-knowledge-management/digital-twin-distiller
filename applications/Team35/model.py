from digital_twin_distiller import DirichletBoundaryCondition, NeumannBoundaryCondition
from digital_twin_distiller import Material
from digital_twin_distiller import BaseModel
from digital_twin_distiller import ModelDir
from digital_twin_distiller import Agros2D, Femm
from digital_twin_distiller import Agros2DMetadata, FemmMetadata
from digital_twin_distiller import Snapshot
from digital_twin_distiller import ModelPiece
import numpy as np

ModelDir.set_base(__file__)

class Team35(BaseModel):
    """docstring for Team35"""
    def __init__(self, X:list, **kwargs):
        assert len(X) == 10
        self.X = X.copy()

        super(Team35, self).__init__(**kwargs)
        self._init_directories()

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "axisymmetric"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False

        agros_metadata = Agros2DMetadata()
        agros_metadata.file_script_name = self.file_solver_script
        agros_metadata.file_metrics_name = self.file_solution
        agros_metadata.problem_type = "magnetic"
        agros_metadata.coordinate_type = "axisymmetric"
        agros_metadata.analysis_type = "steadystate"
        agros_metadata.unit = 1e-3
        agros_metadata.nb_refinements = 0
        agros_metadata.adaptivity = "hp-adaptivity"
        agros_metadata.polyorder = 2
        agros_metadata.adaptivity_tol = 0.001

        platform_femm = Femm(femm_metadata)
        platform_agros = Agros2D(agros_metadata)

        platform = platform_femm
        self.snapshot = Snapshot(platform)

    def define_materials(self):
        exctitation = Material("J+")
        exctitation.Je = 2e6
        # exctitation.meshsize = 1

        air = Material("air")
        # air.meshsize = 0.7

        control = Material("control")
        # control.meshsize = 0.1

        self.snapshot.add_material(exctitation)
        self.snapshot.add_material(air)
        self.snapshot.add_material(control)

    def define_boundary_conditions(self):
        b1 = DirichletBoundaryCondition(name="a0", field_type="magnetic", magnetic_potential=0.0)
        n0 = NeumannBoundaryCondition(name="n0", field_type="magnetic", surface_current=0.0)

        self.snapshot.add_boundary_condition(b1)
        self.snapshot.add_boundary_condition(n0)

        self.boundary_queue.append((0, 2.5, "a0"))
        self.boundary_queue.append((0, 25, "a0"))
        self.boundary_queue.append((40, 140, "a0"))
        self.boundary_queue.append((140, 40, "a0"))

        p0 = (2.5, 0.0)
        p1 = ((self.X[0] + 5) / 2, 0.0)
        p2 = (self.X[0] + 0.5, 0.0)
        p3 = ((self.X[0] + 1.0 + 140) / 2, 0.0)
        
        self.assign_boundary(p0[0], p0[1], "n0")
        self.assign_boundary(p1[0], p1[1], "n0")
        self.assign_boundary(p2[0], p2[1], "n0")
        self.assign_boundary(p3[0], p3[1], "n0")

    def add_postprocessing(self):
        points = [(0, 0)]
        self.snapshot.add_postprocessing("integration", points, "Energy")

    def build_geometry(self):
        # ...
        self.snapshot.add_geometry(self.geom)


if __name__ == "__main__":
    m = Team35(exportname="dev")
    print(m(cleanup=False))
