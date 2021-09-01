from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.boundaries import NeumannBoundaryCondition
from adze_modeler.material import Material
from adze_modeler.metadata import Agros2DMetadata
from adze_modeler.metadata import FemmMetadata
from adze_modeler.model import BaseModel
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot

import numpy as np


class AsymmetircModel(BaseModel):
    """docstring for SymmetircModel"""

    def __init__(self, X, exportname=None):
        assert len(X) == 20
        self.X = X.copy()

        super().__init__(exportname=exportname)
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
        agros_metadata.adaptivity_tol = 1

        platform_femm = Femm(femm_metadata)
        platform_agros = Agros2D(agros_metadata)

        self.platform = platform_femm
        self.snapshot = Snapshot(self.platform)

    def define_boundary_conditions(self):
        b1 = DirichletBoundaryCondition(name="a0", field_type="magnetic", magnetic_potential=0.0)

        self.snapshot.add_boundary_condition(b1)

        self.boundary_queue.append((0, 0, "a0"))
        self.boundary_queue.append((0, -20, "a0"))
        self.boundary_queue.append((0, 20, "a0"))

        self.boundary_queue.append((70, 70, "a0"))
        self.boundary_queue.append((70, -70, "a0"))
        self.boundary_queue.append((140, 0, "a0"))

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

        self.label_queue.append((3, 0, "control"))
        self.label_queue.append((30, 30, "air"))

    def add_postprocessing(self):
        Nx = 10
        Ny = 10
        px = np.linspace(0.001, 5, Nx)
        py = np.linspace(-5, 5, Ny)
        xv, yv = np.meshgrid(px, py, sparse=False, indexing="xy")

        for i in range(Nx):
            for j in range(Ny):
                eval_point = (xv[j, i], yv[j, i])
                self.snapshot.add_postprocessing("point_value", eval_point, "Bz")
                self.snapshot.add_postprocessing("point_value", eval_point, "Br")

        self.snapshot.add_postprocessing("mesh_info", None, None)

    def build_geometry(self):
        mp_bound = ModelPiece("bound")
        mp_bound.load_piece_from_svg(self.dir_resources / "problem_boundary.svg")
        mp_bound.put(0, -70)

        mp_control = ModelPiece("core")
        mp_control.load_piece_from_svg(self.dir_resources / "core_all.svg")
        mp_control.put(0, -5)

        mp_coil = ModelPiece("coil")
        mp_coil.load_piece_from_svg(self.dir_resources / "coil.svg")

        N = len(self.X)
        h = 1.5
        w = 1.0
        self.geom.merge_geometry(mp_bound.geom)
        self.geom.merge_geometry(mp_control.geom)

        for i, ri in enumerate(self.X):
            coil = mp_coil.spawn()
            offsetz = i * h - N * h / 2
            coil.put(ri, offsetz)
            self.geom.merge_geometry(coil.geom)
            self.label_queue.append((ri + w / 2, offsetz + h / 2, "J+"))

        self.geom.generate_intersections()
        self.snapshot.add_geometry(self.geom)


if __name__ == "__main__":
    X = [10] * 20
    m = AsymmetircModel(X)
    print(m(cleanup=False, devmode=False))
