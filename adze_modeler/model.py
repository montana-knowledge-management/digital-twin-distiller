from abc import ABCMeta
from abc import abstractmethod
from types import SimpleNamespace

from adze_modeler.geometry import Geometry


class Snapshot(metaclass=ABCMeta):
    """
    This class contains the full model = geometry + boundary conditions + material parameters.
    The contained data is enough only for a single numerical computation.

    A snapshot can be exported into a
    """

    def __init__(self):
        self.geometry = Geometry
        self.material_definitions = []
        self.labels = []
        self.boundary_conditions = []
        self.epsilon = 1.0e-5

    @abstractmethod
    def export_snapshot(self, fem_code):
        pass


#
# class PhysicalModel:
#     materials = []
#     boundary_conditions = []
#
#
# class Label(SimpleNamespace):
#     def __init__(self, node_id, label_id, material):
#         self.node_id = node_id
#         self.label_id = label_id
#         self.material = material
#
#
# class Material(SimpleNamespace):
#     pass
#
#
# class FemObjects(SimpleNamespace):
#     @abstractmethod
#     def write_geometry(self):
#         pass
#
#
# class Rectangle(FemObjects):
#     def __init__(self, bottom_left: tuple, width, height, material):
#         self.width = width
#         self.height = height
#         self.material = material
#
#
# class BoundaryCondition(SimpleNamespace):
#     def __init__(self, name, type, value):
#         self.name = name
#         self.type = type
#         self.value = value
#

if __name__ == "__main__":
    pass
