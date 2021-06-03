from abc import ABCMeta
from abc import abstractmethod
from types import SimpleNamespace

from adze_modeler.geometry import Geometry


class Snapshot(metaclass=ABCMeta):
    """
    This class contains the full model = geometry + boundary conditions + material parameters.
    The contained data is enough only for a single numerical computation.

    A snapshot can be exported and solved by the defined numerical solver -- as it is --.
    """

    def __init__(self):
        self.geometry = Geometry
        self.region_labels = []
        self.material_definitions = []
        self.boundary_conditions = []
        self.epsilon = 1.0e-5

    def export_snapshot(self, solver):
        """Writes out a single snapshot of the numerical model into the selected FEM library."""
        numerical_solver = solver()
        numerical_solver.write()

    def run(self):
        pass
