from enum import Enum

from digital_twin_distiller.boundaries import DirichletBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.objects import Line, Node
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot

ModelDir.set_base(__file__)


class RatioModel:
    def __init__(self, u1, u2, u3, l1, l2, l3):
        self.u1 = u1
        self.u2 = u2
        self.u3 = u3

        self.l1 = l1
        self.l2 = l2
        self.l3 = l3


class InitialRatioModel:
    def __init__(self, H, L, lambda1, lambda2):
        self.H = H
        self.L = L
        self.lambda1 = lambda1
        self.lambda2 = lambda2


class Boundary(Enum):
    upper = 2500
    grounded = 0


# class SmallRatioModel(Enum):
#     u1 = Node(0, 0.2032)
#     u2 = Node(0, 0.61976)
#     u3 = Node(0.14224, 0.61976)
#
#     l1 = Node(0, 0)
#     l2 = Node(0.14224, 0)
#     l3 = Node(0.14224, 0.41656)
#
#
# class LargeRatioModel(Enum):
#     u1 = Node(0, 0.2032)
#     u2 = Node(0, 1.1176)
#     u3 = Node(0.14224, 1.1176)
#
#     l1 = Node(0, 0)
#     l2 = Node(0.14224, 0)
#     l3 = Node(0.14224, 0.9144)


class SmallModel(InitialRatioModel):

    def __init__(self):
        self.H = 0.07112
        self.L = 0.10668
        self.lambda1 = 0.2032
        self.lambda2 = 0.2032
        super().__init__(self.H, self.L, self.lambda1, self.lambda2)


class LargeModel(InitialRatioModel):

    def __init__(self):
        self.H = 0.07112
        self.L = 0.3556
        self.lambda1 = 0.2032
        self.lambda2 = 0.2032
        super().__init__(self.H, self.L, self.lambda1, self.lambda2)


class SimulationModel(BaseModel):
    """this SimulationModel created when calling 'new'"""

    def __init__(self, **kwargs):
        super(SimulationModel, self).__init__(**kwargs)
        self._init_directories()

    def setup_solver(self):  # TODO: example with FEM solver
        self.platform = Femm(self.define_femm_metadat())
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
        # TODO: create
        ...

    def define_materials(self):
        air = Material('air')
        self.snapshot.add_material(air)

    def define_boundary_conditions(self):
        upper = DirichletBoundaryCondition(Boundary.upper.name, field_type="electrostatic",
                                           fixed_voltage=Boundary.upper.value)
        grounded = DirichletBoundaryCondition(Boundary.grounded.name, field_type="electrostatic",
                                              fixed_voltage=Boundary.grounded.value)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(upper)
        self.snapshot.add_boundary_condition(grounded)

    def add_postprocessing(self):  # TODO
        points = [(0, 0)]
        self.snapshot.add_postprocessing("integration", points, "Energy")

    def build_geometry(self):
        # u1 = Node(0, 0.2032)
        # u2 = Node(0, 0.61976)
        # u3 = Node(0.14224, 0.61976)
        #
        # l1 = Node(0, 0)
        # l2 = Node(0.14224, 0)
        # l3 = Node(0.14224, 0.41656)

        model = self.generate_nodes_with_specific_values(H=0.07112, L=0.10668, lambda1=0.2032, lambda2=0.2032)
        #model = self.generate_nodes_with_predefined_values(model=LargeModel())


        self.build_lines_and_boundaries(model.u1,
                                        model.u2,
                                        model.u3,
                                        model.l1,
                                        model.l2,
                                        model.l3)

        self.snapshot.add_geometry(self.geom)
        self.assign_material(0.1, 0.1, "air")

    # TODO: some coordinates are wrong!
    def generate_nodes_with_specific_values(self, H: float, L: float, lambda1: float, lambda2: float):
        """
        lambda1 and lambda2 are equal, it should be enough to use only one lambda's value
        """
        l1 = Node(0, 0)
        l2 = Node(0, 2 * H)
        l3 = Node(2 * L + lambda1, 2 * H)

        u1 = Node(lambda2, 0)
        u2 = Node(lambda2 + 2 * L + lambda1, 0)
        u3 = Node(2 * H, lambda2 + 2 * L + lambda1)

        return RatioModel(u1=u1, u2=u2, u3=u3, l1=l1, l2=l2, l3=l3)

    def generate_nodes_with_predefined_values(self, model: InitialRatioModel):
        return self.generate_nodes_with_specific_values(H=model.H, L=model.L, lambda1=model.lambda1,
                                                        lambda2=model.lambda2)

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
        grounded = Boundary.grounded.name
        upper = Boundary.upper.name
        self.snapshot.boundaries.get(grounded).assigned.add(line_alpha2.id)
        self.snapshot.boundaries.get(grounded).assigned.add(line_alpha1.id)
        self.snapshot.boundaries.get(grounded).assigned.add(line_l1.id)
        self.snapshot.boundaries.get(grounded).assigned.add(line_l2.id)
        self.snapshot.boundaries.get(upper).assigned.add(line_u1.id)
        self.snapshot.boundaries.get(upper).assigned.add(line_u2.id)


if __name__ == "__main__":
    m = SimulationModel(exportname="dev")
    print(m(cleanup=False, devmode=True))
