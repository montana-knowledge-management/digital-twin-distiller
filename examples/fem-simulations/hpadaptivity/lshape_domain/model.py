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

ModelDir.set_base(__file__)


class LShapeDomain(BaseModel):
    """docstring for lshape_domain"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._init_directories()

        self.l = kwargs.get('l ', 0.001)

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

        self.platform = Agros2D(agros_metadata)

        self.snapshot = Snapshot(self.platform)


    def define_materials(self):
        air = Material("air")
        air.meshsize=0.01

        self.snapshot.add_material(air)

    def define_boundary_conditions(self):
        vcc = DirichletBoundaryCondition("vcc", field_type="electrostatic", fixed_voltage=100.0)
        gnd = DirichletBoundaryCondition("gnd", field_type="electrostatic", fixed_voltage=0.0)

        n0 = NeumannBoundaryCondition("n0", field_type="electrostatic", surface_charge_density=0.0)
        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(vcc)
        self.snapshot.add_boundary_condition(gnd)
        self.snapshot.add_boundary_condition(n0)

    def add_postprocessing(self):
        self.snapshot.add_postprocessing("mesh_info", None, None)

    def build_geometry(self):
        
        p1 = Node(0, -1)
        p2 = Node(1, -1)
        p3 = Node(1, 1)
        p4 = Node(-1, 1)
        p5 = Node(-1, 0)
        p6 = Node(-self.l, 0)
        p7 = Node(0, 0)
        p8 = Node(0, -self.l)

        l1 = Line(p1, p2)
        l2 = Line(p2, p3)
        l3 = Line(p3, p4)
        l4 = Line(p4, p5)
        l5 = Line(p5, p6)
        l6 = Line(p6, p7)
        l7 = Line(p7, p8)
        l8 = Line(p8, p1)

        self.geom.add_line(l1)
        self.geom.add_line(l2)
        self.geom.add_line(l3)
        self.geom.add_line(l4)
        self.geom.add_line(l5)
        self.geom.add_line(l6)
        self.geom.add_line(l7)
        self.geom.add_line(l8)

        self.assign_boundary(*l1(0.5), "n0")
        self.assign_boundary(*l2(0.5), "n0")
        self.assign_boundary(*l3(0.5), "n0")
        self.assign_boundary(*l4(0.5), "n0")
        self.assign_boundary(*l5(0.5), "gnd")
        self.assign_boundary(*l6(0.5), "vcc")
        self.assign_boundary(*l7(0.5), "vcc")
        self.assign_boundary(*l8(0.5), "gnd")

        self.assign_material(0.3, 0.3, "air")

        self.snapshot.add_geometry(self.geom)


if __name__ == "__main__":
    m = LShapeDomain(exportname="dev")
    print(m(cleanup=False, devmode=False))
