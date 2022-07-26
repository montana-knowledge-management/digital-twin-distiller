from copy import deepcopy
from itertools import count
from math import atan, hypot
from os import devnull

from numpy import linspace, meshgrid
from numpy.core.function_base import linspace

from digital_twin_distiller import Agros2D, Agros2DMetadata
from digital_twin_distiller.boundaries import DirichletBoundaryCondition, NeumannBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.objects import Line, Node
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot
from digital_twin_distiller.utils import pairwise

ModelDir.set_base(__file__)

alpha = 50
x0 = 0.5
y0 = 0.5
r0 = 0.25


def u(x, y):
    return atan(alpha * (hypot(x - x0, y - y0) - r0))


def gradu(x, y):
    w = hypot(x - x0, y - y0)
    denominator = w * (1 + 3600 * (w - 1) ** 2)
    dudx = 60 * (x - 1.25) / denominator
    dudy = 60 * (y + 0.25) / denominator
    return dudx, dudy


def gamma(x, y, nx, ny):
    dx, dy = gradu(x, y)
    return dx * nx + dy * ny


class ShockProbem(BaseModel):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._init_directories()

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "electrostatic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "meters"
        femm_metadata.smartmesh = True
        femm_metadata.depth = 1000

        agros_metadata = Agros2DMetadata()
        agros_metadata.file_script_name = self.file_solver_script
        agros_metadata.file_metrics_name = self.file_solution
        agros_metadata.problem_type = "electrostatic"
        agros_metadata.coordinate_type = "planar"
        agros_metadata.analysis_type = "steadystate"
        agros_metadata.unit = 1
        agros_metadata.nb_refinements = 0
        agros_metadata.adaptivity = "hp-adaptivity"
        agros_metadata.polyorder = 4
        agros_metadata.adaptivity_tol = 1
        agros_metadata.adaptivity_steps = 5

        self.platform = Agros2D(agros_metadata)
        self.platform = Femm(femm_metadata)

        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        air = Material("air")
        air.meshsize = 0.1
        air.epsioln_r = 1 / 8.8541878128e-12

        self.snapshot.add_material(air)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="electrostatic", fixed_voltage=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

    def add_postprocessing(self):
        N = 100

        x = linspace(0.5, 1, N)
        y = linspace(0.5, 1, N)
        xx, yy = meshgrid(x, y)

        for i in range(N):
            for j in range(N):
                eval_point = (xx[j, i], yy[j, i])
                self.snapshot.add_postprocessing("point_value", eval_point, "V")

    def build_geometry(self):
        N = 500
        r0 = Node(0.5, 0.5)
        r1 = Node(1, 0.5)
        r2 = Node(1, 1)
        r3 = Node(0.5, 1)

        # g1 = lambda x, y: gamma(x, y, 0, 1)
        # g2 = lambda x, y: gamma(x, y, 1, 0)
        # g3 = lambda x, y: gamma(x, y, 0, 1)
        # g4 = lambda x, y: gamma(x, y, -1, 0)

        n0 = NeumannBoundaryCondition("gamma_1", field_type="electrostatic")
        n1 = NeumannBoundaryCondition("gamma_2", field_type="electrostatic")
        n2 = NeumannBoundaryCondition("gamma_3", field_type="electrostatic")
        n3 = NeumannBoundaryCondition("gamma_4", field_type="electrostatic")

        a0 = DirichletBoundaryCondition("gamma_D_1", field_type="electrostatic")
        a1 = DirichletBoundaryCondition("gamma_D_2", field_type="electrostatic")
        a2 = DirichletBoundaryCondition("gamma_D_3", field_type="electrostatic")
        a3 = DirichletBoundaryCondition("gamma_D_4", field_type="electrostatic")

        self._approximate_boundary(a0, "fixed_voltage", r0, r1, N, u)
        self._approximate_boundary(a1, "fixed_voltage", r1, r2, N, u)
        self._approximate_boundary(a2, "fixed_voltage", r2, r3, N, u)
        self._approximate_boundary(a3, "fixed_voltage", r3, r0, N, u)

        self.assign_material(0.95, 0.95, "air")

        self.snapshot.add_geometry(self.geom)

    def _approximate_boundary(self, b0, variable, p0, p1, N, fun_):
        l = Line(p0, p1)
        t = linspace(0, 1, N + 1)
        index = count(start=1)
        for ti, tj in pairwise(t):
            bi = deepcopy(b0)
            bi.name = f"{bi.name}_{next(index)}"
            segment = Line(l(ti), l(tj))
            center_point = segment(0.5)
            boundary_value = fun_(center_point.x, center_point.y)
            bi.set_value(variable, boundary_value)
            bi.assigned.add(segment.id)
            self.snapshot.add_boundary_condition(bi)

            self.geom.nodes.append(segment.start_pt)
            self.geom.nodes.append(segment.end_pt)
            self.geom.lines.append(segment)


if __name__ == "__main__":
    import csv

    m = ShockProbem(exportname="dev")
    r_ = m(cleanup=False, devmode=False)

    with open(ModelDir.DATA / "r_femm.csv", "w") as f:
        w = csv.writer(f)
        for ri in r_["V"]:
            w.writerow(ri)
