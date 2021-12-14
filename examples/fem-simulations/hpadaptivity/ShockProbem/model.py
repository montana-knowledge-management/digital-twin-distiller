from itertools import count
from copy import deepcopy
from math import atan, hypot
from os import devnull
from numpy.core.function_base import linspace
from digital_twin_distiller.boundaries import DirichletBoundaryCondition, NeumannBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.objects import Line, Node
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot
from digital_twin_distiller import Agros2DMetadata, Agros2D
from digital_twin_distiller.utils import pairwise

ModelDir.set_base(__file__)

def u(x, y):
    return atan(60*(hypot(x-1.25, y+0.25)-1))

def gradu(x, y):
    w = hypot(x-1.25, y+0.25)
    denominator = w * (1+3600*(w-1)**2)
    dudx = 60 * (x - 1.25) / denominator
    dudy = 60 * (y + 0.25) / denominator
    return dudx, dudy

def gamma(x, y, nx, ny):
    dx, dy = gradu(x, y)
    return dx * nx + dy * ny

class ShockProbem(BaseModel):
    def __init__(self, **kwargs):
        super(ShockProbem, self).__init__(**kwargs)
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
        agros_metadata.adaptivity_tol = 10
        agros_metadata.adaptivity_steps = 50

        self.platform = Agros2D(agros_metadata)
        self.platform = Femm(femm_metadata)


        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        air = Material("air")
        air.meshsize = 0.1
        air.epsioln_r = (1/8.8541878128e-12)

        self.snapshot.add_material(air)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="electrostatic", fixed_voltage=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

    def add_postprocessing(self):
        points = [(0, 0)]
        self.snapshot.add_postprocessing("integration", points, "Energy")

    def build_geometry(self):
        N = 200
        r0 = Node(0, 0)
        r1 = Node(1, 0)
        r2 = Node(1, 1)
        r3 = Node(0, 1)

        g1 = lambda x, y: gamma(x, y, 0, 1)
        g2 = lambda x, y: gamma(x, y, 1, 0)
        g3 = lambda x, y: gamma(x, y, 0, 1)
        g4 = lambda x, y: gamma(x, y, -1, 0)

        n0 = NeumannBoundaryCondition("gamma_1", field_type="electrostatic")
        n1 = NeumannBoundaryCondition("gamma_2", field_type="electrostatic")
        n2 = NeumannBoundaryCondition("gamma_3", field_type="electrostatic")
        n3 = NeumannBoundaryCondition("gamma_4", field_type="electrostatic")

        a0 = DirichletBoundaryCondition("gamma_D_1", field_type="electrostatic")
        a1 = DirichletBoundaryCondition("gamma_D_2", field_type="electrostatic")
        a2 = DirichletBoundaryCondition("gamma_D_3", field_type="electrostatic")
        a3 = DirichletBoundaryCondition("gamma_D_4", field_type="electrostatic")

        self._approximate_boundary(n0, "surface_charge_density", r0, r1, N, g1)
        self._approximate_boundary(n1, "surface_charge_density", r1, r2, N, g2)
        # self._approximate_boundary(n2, "surface_charge_density", r2, r3, N, g3)
        # self._approximate_boundary(n3, "surface_charge_density", r3, r0, N, g4)
        self._approximate_boundary(a2, "fixed_voltage", r2, r3, N, u)
        self._approximate_boundary(a3, "fixed_voltage", r3, r0, N, u)


        self.assign_material(0.5, 0.5, "air")

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
    m = ShockProbem(exportname="dev")
    print(m(cleanup=False, devmode=False))
