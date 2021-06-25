from abc import ABCMeta
from abc import abstractmethod
from adze_modeler.material import Material
from adze_modeler.boundaries import BoundaryCondition
from adze_modeler.objects import Node, Line, CircleArc, CubicBezier
from adze_modeler.metadata import Metadata

class Platform(metaclass=ABCMeta):
    def __init__(self, m: Metadata):
        self.metadata = m
        self.metadata.validate_metadata()

        self.file_script_handle = None



    def open(self):
        self.file_script_handle = open(self.metadata.file_script_name, "w")

    def close(self):
        self.file_script_handle.close()

    def write(self, str_, nb_newline=1):
        self.file_script_handle.write(str_)
        self.newline(nb_newline)

    def newline(self, n):
        self.file_script_handle.write('\n' * n)


    def get_handle(self):
        return self.file_script_handle

    def get_script_name(self):
        filename = self.metadata.file_script_name
        dotindex = filename.find('.')
        filename = filename[:dotindex]
        return filename

    @abstractmethod
    def comment(self, str_, nb_newline=1):
        ...

    @abstractmethod
    def export_preamble(self):
        ...

    @abstractmethod
    def export_metadata(self):
        ...

    @abstractmethod
    def export_material_definition(self, mat: Material):
        ...

    @abstractmethod
    def export_block_label(self, x, y, mat: Material):
        ...

    @abstractmethod
    def export_boundary_definition(self, b: BoundaryCondition):
        ...

    @abstractmethod
    def export_geometry_element(self, e, boundary=None):
        ...

    @abstractmethod
    def export_solving_steps(self):
        ...


    @abstractmethod
    def export_metrics(self, action, entity, variable):
        ...

    @abstractmethod
    def export_closing_steps(self):
        ...

    @abstractmethod
    def execute(self):
        ...