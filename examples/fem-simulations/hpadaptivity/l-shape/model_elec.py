from itertools import count
from typing import Callable
from numpy import linspace
from digital_twin_distiller.boundaries import DirichletBoundaryCondition
from digital_twin_distiller import NeumannBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.objects import Line, Node
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller import Agros2D, Agros2DMetadata
from digital_twin_distiller.snapshot import Snapshot
from math import atan, cos, hypot, atan2, sin, pi, sqrt

from digital_twin_distiller.utils import pairwise

ModelDir.set_base(__file__)

X, Y = 0.7, 0.0


def u(x, y):
    r = hypot(x, y)
    phi = atan2(y, x)
    return r ** (2 / 3) * sin(2 / 3 * phi + pi / 3)


def gradu(x, y):
    alpha = 2 * atan2(y, x) / 3 + pi / 3
    dudx = (2 * sin(alpha) * x - 2 * cos(alpha) * y) / (3 * hypot(x, y) ** (2 / 3))
    dudy = (2 * sin(alpha) * y + 2 * cos(alpha) * x) / (3 * hypot(x, y) ** (2 / 3))
    return (dudx, dudy)


def gamma(x, y, xu, yu):
    dudx, dudy = gradu(x, y)
    return dudx * xu + dudy * yu

class LShape(BaseModel):
    def __init__(self, **kwargs):
        super(LShape, self).__init__(**kwargs)
        self._init_directories()

    def define_materials(self):
        air = Material("air")
        air.meshsize = 0.01
        air.epsioln_r = (1/8.856e-12)

        self.snapshot.add_material(air)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="electrostatic", fixed_voltage=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

    def add_postprocessing(self):
        # self.snapshot.add_postprocessing("point_value", (0, 0), "T")
        self.snapshot.add_postprocessing("point_value", (X, Y), "V")
        self.snapshot.add_postprocessing("mesh_info", None, None)


class FemmLShape(LShape):
    def __init__(self, **kwargs):
        super(FemmLShape, self).__init__(**kwargs)

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "electrostatic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "meters"
        femm_metadata.smartmesh = False
        femm_metadata.depth = 1

        self.platform = Femm(femm_metadata)

        self.snapshot = Snapshot(self.platform)

    def define_boundary_conditions(self):
        super(FemmLShape, self).define_boundary_conditions()

    def build_geometry(self):
        N = 200

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

        self._approximate_boundary("gamma_1", "fixed_voltage", r1, r2, 2 * N, lambda x, y: gamma(x, y, 0, 1))
        self._approximate_boundary("gamma_2", "fixed_voltage", r2, r3, 2 * N, lambda x, y: gamma(x, y, 1, 0))
        self._approximate_boundary("gamma_3", "fixed_voltage", r3, r4, N, lambda x, y: gamma(x, y, 0, -1))
        self._approximate_boundary("gamma_6", "fixed_voltage", r6, r1, N, lambda x, y: gamma(x, y, -1, 0))

        self.assign_material(0.3, 0.3, "air")

        self.snapshot.add_geometry(self.geom)

    def _approximate_boundary(
        self,
        boundary_prefix: str,
        boundary_variable: str,
        p1: Node,
        p2: Node,
        N: int,
        fun_: Callable[[float, float], float],
    ):
        l = Line(p1, p2)
        t = linspace(0, 1, N + 1)
        index = count(start=1)
        for ti, tj in pairwise(t):
            boundary_id = f"{boundary_prefix}_{next(index)}"
            segment = Line(l(ti), l(tj))
            center_point = segment(0.5)
            boundary_value = fun_(center_point.x, center_point.y)

            ni = NeumannBoundaryCondition(boundary_id, field_type="electrostatic", surface_charge_density=boundary_value)
            ni.assigned.add(segment.id)
            self.snapshot.add_boundary_condition(ni)

            self.geom.nodes.append(segment.start_pt)
            self.geom.nodes.append(segment.end_pt)
            self.geom.lines.append(segment)


class AgrosLShape(LShape):
    def __init__(self, **kwargs):
        super(AgrosLShape, self).__init__(**kwargs)

    def setup_solver(self):
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
        agros_metadata.adaptivity_tol = 0.001
        agros_metadata.adaptivity_steps = 50

        self.platform = Agros2D(agros_metadata)

        self.snapshot = Snapshot(self.platform)

    def define_boundary_conditions(self):
        super(AgrosLShape, self).define_boundary_conditions()

    def build_geometry(self):
        N = 200

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

        self._approximate_boundary("gamma_1", "fixed_voltage", r1, r2, 2 * N, lambda x, y: gamma(x, y, 0, 1))
        self._approximate_boundary("gamma_2", "fixed_voltage", r2, r3, 2 * N, lambda x, y: gamma(x, y, 1, 0))
        self._approximate_boundary("gamma_3", "fixed_voltage", r3, r4, N, lambda x, y: gamma(x, y, 0, -1))
        self._approximate_boundary("gamma_6", "fixed_voltage", r6, r1, N, lambda x, y: gamma(x, y, -1, 0))

        self.assign_material(0.3, 0.3, "air")

        self.snapshot.add_geometry(self.geom)

    def _approximate_boundary(
        self,
        boundary_prefix: str,
        boundary_variable: str,
        p1: Node,
        p2: Node,
        N: int,
        fun_: Callable[[float, float], float],
    ):
        l = Line(p1, p2)
        t = linspace(0, 1, N + 1)
        index = count(start=1)
        for ti, tj in pairwise(t):
            boundary_id = f"{boundary_prefix}_{next(index)}"
            segment = Line(l(ti), l(tj))
            center_point = segment(0.5)
            boundary_value = fun_(center_point.x, center_point.y)

            ni = NeumannBoundaryCondition(boundary_id, field_type="electrostatic", surface_charge_density=boundary_value)
            ni.assigned.add(segment.id)
            self.snapshot.add_boundary_condition(ni)

            self.geom.nodes.append(segment.start_pt)
            self.geom.nodes.append(segment.end_pt)
            self.geom.lines.append(segment)

if __name__ == "__main__":
    m = AgrosLShape(exportname="dev_agros")
    agres = m(cleanup=False, devmode=False)
    print(*agres["V"])

    m = FemmLShape(exportname="dev_femm")
    femres = m(cleanup=False, devmode=False)
    print(*femres["V"])

    print(u(X, Y))
    """
    dudx: (2*sin(2/3*atan2(y,x)+pi/3)*x-2*cos(2/3*atan2(y,x)+pi/3)*y)/(3*hypot(x,y)**(2/3))
    dudy: (2*sin(2/3*atan2(y,x)+pi/3)*y+2*cos(2/3*atan2(y,x)+pi/3)*x)/(3*hypot(x,y)**(2/3))
    """
