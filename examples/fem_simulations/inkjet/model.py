from enum import Enum

from digital_twin_distiller import Agros2D, Platform
from digital_twin_distiller.boundaries import DirichletBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata, Agros2DMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.objects import Line, Node
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot

ModelDir.set_base(__file__)


class NodeModel:
    def __init__(self, u1: Node, u2: Node, u3: Node, l1: Node, l2: Node, l3: Node):
        self.u1 = u1
        self.u2 = u2
        self.u3 = u3

        self.l1 = l1
        self.l2 = l2
        self.l3 = l3


class ParamModel:
    node: NodeModel = None

    def __init__(self, H: float, L: float, lambda1: float, lambda2: float):
        self.H = H
        self.L = L
        self.lambda1 = lambda1
        self.lambda2 = lambda2
        self.init_nodes()

    def init_nodes(self):
        return self.generate_nodes_with_specific_values(H=self.H, L=self.L,
                                                        lambda1=self.lambda1,
                                                        lambda2=self.lambda2)

    def generate_nodes_with_specific_values(self, H: float, L: float, lambda1: float, lambda2: float):
        """
        lambda1 and lambda2 are equal, it should be enough to use only one lambda's value
        """
        l1 = Node(0, 0)
        l2 = Node(2 * H, 0)
        l3 = Node(2 * H, 2 * L + lambda1)

        u1 = Node(0, lambda2)
        u2 = Node(0, lambda2 + 2 * L + lambda1)
        u3 = Node(2 * H, lambda2 + 2 * L + lambda1)

        self.node = NodeModel(u1=u1, u2=u2, u3=u3, l1=l1, l2=l2, l3=l3)


def create_model(H: float, L: float, lambda1: float, lambda2: float):
    return ParamModel(H=H, L=L, lambda1=lambda1, lambda2=lambda2)


def large_model():
    H = 0.07112
    L = 0.3556
    lambda1 = 0.2032
    lambda2 = 0.2032
    return create_model(H=H, L=L, lambda1=lambda1, lambda2=lambda2)


def small_model():
    H = 0.07112
    L = 0.10668
    lambda1 = 0.2032
    lambda2 = 0.2032
    return create_model(H=H, L=L, lambda1=lambda1, lambda2=lambda2)


class PlatformVariable:
    AGROS2D = "agros2d"
    FEMM = "femm"


class SimulationModel(BaseModel):
    """this SimulationModel created when calling 'new'"""

    class Boundary(Enum):
        upper = 2500
        ground = 0

    model: ParamModel = None
    # solver: Platform = None
    ratio: float = None

    def __init__(self, model_params=None, solver=None, **kwargs):
        super(SimulationModel, self).__init__(**kwargs)

        if model_params is None:
            self.model = large_model()
        else:
            self.model = model_params
        if solver == PlatformVariable.AGROS2D:
            self.platform = Agros2D(self.define_agros2d_metadata())
        elif solver == PlatformVariable.FEMM or solver is None:
            self.platform = Femm(self.define_femm_metadat())
        self.calculate_ratio()
        print("Ratio(L/H): ", self.ratio)
        self._init_directories()

    def setup_solver(self):
        # platform_agros = Agros2D(self.define_agros2d_metadata())
        # platform_femm = Femm(self.define_femm_metadat())
        # self.platform = self.solver
        self.snapshot = Snapshot(self.platform)

    def define_femm_metadat(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "electrostatic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "centimeters"
        femm_metadata.smartmesh = True
        femm_metadata.depth = 1000
        return femm_metadata

    def define_agros2d_metadata(self):
        agros_metadata = Agros2DMetadata()
        agros_metadata.file_script_name = self.file_solver_script
        agros_metadata.file_metrics_name = self.file_solution
        agros_metadata.problem_type = "electrostatic"
        agros_metadata.coordinate_type = "planar"
        agros_metadata.analysis_type = "steadystate"
        agros_metadata.unit = 1e-3
        agros_metadata.nb_refinements = 2  # 4
        agros_metadata.adaptivity = "hp-adaptivity"
        agros_metadata.polyorder = 5
        agros_metadata.adaptivity_tol = 1
        return agros_metadata

    def define_materials(self):
        air = Material('air')
        self.snapshot.add_material(air)

    def define_boundary_conditions(self):
        upper = DirichletBoundaryCondition(self.Boundary.upper.name, field_type="electrostatic",
                                           fixed_voltage=self.Boundary.upper.value)
        ground = DirichletBoundaryCondition(self.Boundary.ground.name, field_type="electrostatic",
                                            fixed_voltage=self.Boundary.ground.value)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(upper)
        self.snapshot.add_boundary_condition(ground)

    def add_postprocessing(self):  # TODO
        points = [(0, 0)]
        self.snapshot.add_postprocessing("integration", points, "Energy")

    def build_geometry(self):

        self.build_lines_and_boundaries(self.model.node.u1,
                                        self.model.node.u2,
                                        self.model.node.u3,
                                        self.model.node.l1,
                                        self.model.node.l2,
                                        self.model.node.l3)

        self.snapshot.add_geometry(self.geom)
        self.assign_material(0.1, 0.1, "air")

    def add_geom_lines(self, lines):
        for line in lines:
            self.geom.add_line(line)

    def build_lines_and_boundaries(self, u1: Node, u2: Node, u3: Node, l1: Node, l2: Node, l3: Node):
        line_u1 = Line(u1, u2)
        line_u2 = Line(u2, u3)

        line_l1 = Line(l1, l2)
        line_l2 = Line(l2, l3)

        line_alpha2 = Line(l1, u1)
        line_alpha1 = Line(l3, u3)

        lines = [line_u1, line_u2, line_l1, line_l2, line_alpha2, line_alpha1]
        self.add_geom_lines(lines)

        # connect boundary conditions to lines
        ground = self.Boundary.ground.name
        upper = self.Boundary.upper.name
        self.snapshot.boundaries.get(ground).assigned.add(line_alpha2.id)
        self.snapshot.boundaries.get(ground).assigned.add(line_alpha1.id)
        self.snapshot.boundaries.get(ground).assigned.add(line_l1.id)
        self.snapshot.boundaries.get(ground).assigned.add(line_l2.id)
        self.snapshot.boundaries.get(upper).assigned.add(line_u1.id)
        self.snapshot.boundaries.get(upper).assigned.add(line_u2.id)

    def calculate_ratio(self):
        if self.model.H != 0:
            self.ratio = self.model.L / self.model.H
        else:
            print("H cannot be 0")
            exit(1)


# TODO: create tests
if __name__ == "__main__":
    platform = PlatformVariable.AGROS2D
    # model = create_model(H=0.07112, L=0, lambda1=0.2032, lambda2=0.2032)
    model = small_model()
    # model = large_model()  # default if model param is not given

    m = SimulationModel(model, platform, exportname="dev")
    print(m(cleanup=False, devmode=False))
