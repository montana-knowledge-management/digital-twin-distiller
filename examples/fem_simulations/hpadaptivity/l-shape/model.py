from copy import deepcopy
from itertools import count
from math import atan, atan2, cos, hypot, pi, prod, sin, sqrt
from typing import Callable

from numpy import linspace

from digital_twin_distiller import Agros2D, Agros2DMetadata, NeumannBoundaryCondition
from digital_twin_distiller.boundaries import DirichletBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.objects import Line, Node
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot
from digital_twin_distiller.utils import pairwise

ModelDir.set_base(__file__)


def u(x, y):
    r = hypot(x, y)
    phi = atan2(y, x)
    return r ** (2 / 3) * sin(2 / 3 * phi + pi / 3)


def gradu(x, y):
    alpha = 2 * atan2(y, x) / 3 + pi / 3
    dudx = (2 * sin(alpha) * x - 2 * cos(alpha) * y) / (3 * (x**2 + y**2) ** (2 / 3))
    dudy = (2 * sin(alpha) * y + 2 * cos(alpha) * x) / (3 * (x**2 + y**2) ** (2 / 3))
    return (dudx, dudy)


def gamma(x, y, xu, yu):
    dudx, dudy = gradu(x, y)
    return dudx * xu + dudy * yu


class LShape(BaseModel):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._init_directories()

        self.solver = kwargs.get("solver", "femm")
        self.nb_segments = kwargs.get("nb_segments", 100)

        self.meshsize = kwargs.get("meshsize", 1.345)
        self.adaptivity_tol = kwargs.get("adaptivity_tol", 1)

    def setup_solver(self):
        if self.solver == "femm":
            femm_metadata = FemmMetadata()
            femm_metadata.problem_type = "electrostatic"
            femm_metadata.coordinate_type = "planar"
            femm_metadata.file_script_name = self.file_solver_script
            femm_metadata.file_metrics_name = self.file_solution
            femm_metadata.unit = "meters"
            femm_metadata.smartmesh = False
            femm_metadata.depth = 1

            self.platform = Femm(femm_metadata)
        elif self.solver == "agros2d":
            agros_metadata = Agros2DMetadata()
            agros_metadata.file_script_name = self.file_solver_script
            agros_metadata.file_metrics_name = self.file_solution
            agros_metadata.problem_type = "electrostatic"
            agros_metadata.coordinate_type = "planar"
            agros_metadata.analysis_type = "steadystate"
            agros_metadata.unit = 1
            agros_metadata.nb_refinements = 0
            agros_metadata.adaptivity = "hp-adaptivity"
            agros_metadata.polyorder = 1
            agros_metadata.adaptivity_tol = self.adaptivity_tol
            agros_metadata.adaptivity_steps = 99

            self.platform = Agros2D(agros_metadata)

        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        air = Material("air")
        air.meshsize = self.meshsize
        air.epsioln_r = 1 / 8.8541878128e-12

        self.snapshot.add_material(air)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="electrostatic", fixed_voltage=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

    def build_geometry(self):

        r1 = Node(-1, 1)
        r2 = Node(1, 1)
        r3 = Node(1, -1)
        r4 = Node(0, -1)
        r5 = Node(0, 0)
        r6 = Node(-1, 0)

        l4 = Line(r4, r5)
        l5 = Line(r5, r6)

        self.geom.add_line(l4)
        self.geom.add_line(l5)
        self.snapshot.boundaries.get("a0").assigned.add(l4.id)
        self.snapshot.boundaries.get("a0").assigned.add(l5.id)

        n1 = NeumannBoundaryCondition("gamma_1", field_type="electrostatic")
        n2 = NeumannBoundaryCondition("gamma_2", field_type="electrostatic")
        n3 = NeumannBoundaryCondition("gamma_3", field_type="electrostatic")
        n6 = NeumannBoundaryCondition("gamma_6", field_type="electrostatic")

        g1 = lambda x, y: gamma(x, y, 0, 1)
        g2 = lambda x, y: gamma(x, y, 1, 0)
        g3 = lambda x, y: gamma(x, y, 0, -1)
        g6 = lambda x, y: gamma(x, y, -1, 0)

        self._approximate_boundary(n1, "surface_charge_density", r1, r2, self.nb_segments, g1)
        self._approximate_boundary(n2, "surface_charge_density", r2, r3, self.nb_segments, g2)
        self._approximate_boundary(n3, "surface_charge_density", r3, r4, self.nb_segments // 2, g3)
        self._approximate_boundary(n6, "surface_charge_density", r6, r1, self.nb_segments // 2, g6)

        self.assign_material(0.3, 0.3, "air")

        self.snapshot.add_geometry(self.geom)

    def add_postprocessing(self):
        dx = 1e-9
        dy = 1e-9

        self.snapshot.add_postprocessing("point_value", (0.0, 0.0), "V")
        self.snapshot.add_postprocessing("point_value", (dx, dy), "V")
        self.snapshot.add_postprocessing("point_value", (0.95, 0.95), "V")

        self.snapshot.add_postprocessing("point_value", (0, 0), "Ex")
        self.snapshot.add_postprocessing("point_value", (0, 0), "Ey")

        self.snapshot.add_postprocessing("point_value", (dx, dy), "Ex")
        self.snapshot.add_postprocessing("point_value", (dx, dy), "Ey")

        self.snapshot.add_postprocessing("mesh_info", None, None)

        if self.solver == "agros2d":
            self.snapshot.add_postprocessing("integration", [0], "Energy")
        else:
            self.snapshot.add_postprocessing("integration", [(0.1, 0.1)], "Energy")

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

    m = LShape(exportname="dev", meshsize=0.01, solver="femm")
    res = m(cleanup=False, devmode=False)
    import pprint

    pprint.pprint(res)
    print(u(0.95, 0.95))
