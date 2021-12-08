from itertools import count
from math import atan, hypot, sqrt

from numpy import arctan, linspace

from digital_twin_distiller import Line, Node
from digital_twin_distiller.boundaries import DirichletBoundaryCondition, NeumannBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot
from digital_twin_distiller.utils import pairwise
from digital_twin_distiller import Agros2DMetadata, Agros2D

ModelDir.set_base(__file__)


class ThermalShock(BaseModel):
    """docstring for hpadaptivity"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._init_directories()

        self.alpha = kwargs.get("alpha", 60.0)
        self.r0 = kwargs.get("r0", 1.0)

        self.nsteps = kwargs.get("nsteps", 10)

        self.n_gamma1 = kwargs.get("n_gamma1", self.nsteps)
        self.n_gamma2 = kwargs.get("n_gamma2", self.nsteps)
        self.n_gamma3 = kwargs.get("n_gamma3", self.nsteps)
        self.n_gamma4 = kwargs.get("n_gamma4", self.nsteps)

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
        agros_metadata.adaptivity_tol = 0.15
        agros_metadata.adaptivity_steps = 10
        
        # self.platform = Agros2D(agros_metadata)

        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        air = Material("air")
        air.meshsize = 0.005

        self.snapshot.add_material(air)

    def define_boundary_conditions(self):

        for i in range(self.n_gamma1):
            gamma_i = DirichletBoundaryCondition(f"gamma_1_{i+1}", field_type="electrostatic", fixed_voltage=0.0)
            self.snapshot.add_boundary_condition(gamma_i)

        for i in range(self.n_gamma2):
            gamma_i = DirichletBoundaryCondition(f"gamma_2_{i+1}", field_type="electrostatic", fixed_voltage=0.0)
            self.snapshot.add_boundary_condition(gamma_i)

        for i in range(self.n_gamma3):
            gamma_i = DirichletBoundaryCondition(f"gamma_3_{i+1}", field_type="electrostatic", fixed_voltage=0.0)
            self.snapshot.add_boundary_condition(gamma_i)

        for i in range(self.n_gamma4):
            gamma_i = DirichletBoundaryCondition(f"gamma_4_{i+1}", field_type="electrostatic", fixed_voltage=0.0)
            self.snapshot.add_boundary_condition(gamma_i)

        n0 = NeumannBoundaryCondition("n0", field_type="electrostatic", surface_charge_density=0.0)
        self.snapshot.add_boundary_condition(n0)

    def add_postprocessing(self):
        self.snapshot.add_postprocessing("mesh_info", None, None)

    def build_geometry(self):
        # ...
        a = Node(0, 0)
        b = Node(1, 0)
        c = Node(1, 1)
        d = Node(0, 1)

        l1 = Line(a, b)
        l2 = Line(b, c)
        l3 = Line(d, c)
        l4 = Line(a, d)
        self.assign_material(0.5, 0.5, "air")

        u = lambda x, y: arctan(self.alpha*(hypot(1.25-x, -0.25-y)-1))

        # GAMMA 1
        t1 = linspace(0, 1, self.n_gamma1 + 1)
        index = count(start=1)
        for ti, tj in pairwise(t1):
            boundary_id = f"gamma_1_{next(index)}"
            l = Line(l1(ti), l1(tj))
            center_point = l(0.5)
            gamma_value = u(*center_point)
            self.snapshot.boundaries.get(boundary_id).set_value("fixed_voltage", gamma_value)
            self.geom.add_line(l)
            # self.assign_boundary(center_point.x, center_point.y, boundary_id)
            # self.geom.nodes.append(l.start_pt)
            # self.geom.nodes.append(l.end_pt)
            # self.geom.lines.append(l)
            self.snapshot.boundaries.get(boundary_id).assigned.add(l.id)


        # GAMMA 2
        t2 = linspace(0, 1, self.n_gamma2 + 1)
        index = count(start=1)
        for ti, tj in pairwise(t2):
            boundary_id = f"gamma_2_{next(index)}"
            l = Line(l2(ti), l2(tj))
            center_point = l(0.5)
            gamma_value = u(*center_point)
            self.snapshot.boundaries.get(boundary_id).set_value("fixed_voltage", gamma_value)
            self.geom.add_line(l)
            # self.assign_boundary(center_point.x, center_point.y, boundary_id)

            # self.geom.nodes.append(l.start_pt)
            # self.geom.nodes.append(l.end_pt)
            # self.geom.lines.append(l)
            self.snapshot.boundaries.get(boundary_id).assigned.add(l.id)

        # self.geom.nodes.append(l2.start_pt)
        # self.geom.nodes.append(l2.end_pt)
        # self.geom.lines.append(l2)
        # self.snapshot.boundaries.get("n0").assigned.add(l2.id)


        # GAMMA 3
        t3 = linspace(0, 1, self.n_gamma3 + 1)
        index = count(start=1)
        for ti, tj in pairwise(t3):
            boundary_id = f"gamma_3_{next(index)}"
            l = Line(l3(ti), l3(tj))
            center_point = l(0.5)
            gamma_value = u(*center_point)
            self.snapshot.boundaries.get(boundary_id).set_value("fixed_voltage", gamma_value)
            self.geom.add_line(l)
            # self.assign_boundary(center_point.x, center_point.y, boundary_id)

            # self.geom.nodes.append(l.start_pt)
            # self.geom.nodes.append(l.end_pt)
            # self.geom.lines.append(l)
            self.snapshot.boundaries.get(boundary_id).assigned.add(l.id)

        # GAMMA 4
        t4 = linspace(0, 1, self.n_gamma4 + 1)
        index = count(start=1)
        for ti, tj in pairwise(t4):
            boundary_id = f"gamma_4_{next(index)}"
            l = Line(l4(ti), l4(tj))
            center_point = l(0.5)
            gamma_value = u(*center_point)
            self.snapshot.boundaries.get(boundary_id).set_value("fixed_voltage", gamma_value)
            self.geom.add_line(l)
            # self.assign_boundary(center_point.x, center_point.y, boundary_id)

            # self.geom.nodes.append(l.start_pt)
            # self.geom.nodes.append(l.end_pt)
            # self.geom.lines.append(l)
            self.snapshot.boundaries.get(boundary_id).assigned.add(l.id)

        self.snapshot.add_geometry(self.geom)


if __name__ == "__main__":
    m = ThermalShock(exportname="dev", nsteps=100)
    print(m(cleanup=False, devmode=False))
