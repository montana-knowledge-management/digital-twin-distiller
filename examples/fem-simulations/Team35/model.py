import numpy as np

from digital_twin_distiller import (
    Agros2D,
    Agros2DMetadata,
    BaseModel,
    DirichletBoundaryCondition,
    Femm,
    FemmMetadata,
    Material,
    ModelDir,
    ModelPiece,
    NeumannBoundaryCondition,
    Snapshot,
)

ModelDir.set_base(__file__)


class DistributedWinding(BaseModel):
    """docstring for Team35"""

    def __init__(self, X: list, **kwargs):
        assert len(X) == 10
        self.X = X.copy()

        super().__init__(**kwargs)
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
        # agros_metadata.adaptivity_tol = 0.001
        agros_metadata.adaptivity_tol = 1

        platform_femm = Femm(femm_metadata)
        platform_agros = Agros2D(agros_metadata)

        platform = platform_femm
        # platform = platform_agros
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
        Nx = 10
        Ny = 10
        px = np.linspace(0.001, 5, Nx)
        py = np.linspace(0.001, 5, Ny)
        xv, yv = np.meshgrid(px, py, sparse=False, indexing="xy")

        for i in range(Nx):
            for j in range(Ny):
                eval_point = (xv[j, i], yv[j, i])
                self.snapshot.add_postprocessing("point_value", eval_point, "Bz")
                self.snapshot.add_postprocessing("point_value", eval_point, "Br")

        self.snapshot.add_postprocessing("mesh_info", None, None)

    def build_geometry(self):
        mp_bound = ModelPiece("bound")
        mp_bound.load_piece_from_svg(ModelDir.RESOURCES / "problem_boundary.svg")
        mp_bound.put(0, 0)

        mp_control = ModelPiece("core")
        mp_control.load_piece_from_svg(ModelDir.RESOURCES / "core_half.svg")
        mp_control.put(0, 0)

        mp_coil = ModelPiece("coil")
        mp_coil.load_piece_from_svg(ModelDir.RESOURCES / "coil.svg")

        h = 1.5
        w = 1.0
        self.geom.merge_geometry(mp_bound.geom)
        self.geom.merge_geometry(mp_control.geom)

        for i, ri in enumerate(self.X):
            coil = mp_coil.spawn()
            offsetz = i * h
            coil.put(ri, offsetz)
            self.geom.merge_geometry(coil.geom)
            self.label_queue.append((ri + w / 2, offsetz + h / 2, "J+"))

        self.geom.generate_intersections()
        self.snapshot.add_geometry(self.geom)

        self.assign_material(3, 1, "control")
        self.assign_material(30, 30, "air")


if __name__ == "__main__":
    X = [10] * 10
    m = DistributedWinding(X, exportname="dev")
    print(m(cleanup=False, devmode=True))
