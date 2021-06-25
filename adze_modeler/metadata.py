from abc import ABCMeta, abstractmethod


class Metadata(metaclass=ABCMeta):
    def __init__(self):
        self.compatible_platform = None
        self.problem_type = None
        self.analysis_type = None
        self.coordinate_type = None
        self.mesh_type = None
        self.unit = 1.0 # 1 m
        self.precision = 1e-8

        self.file_suffix = None
        self.file_script_name = None
        self.file_metrics_name = "fem_data.csv"

    def validate_file_name(self):
        dotindex = self.file_script_name.find('.')
        if dotindex != -1:
            self.file_script_name = self.file_script_name[:dotindex]

        self.file_script_name = self.file_script_name + self.file_suffix

    @abstractmethod
    def validate_metadata(self):
        ...



class Agros2DMetadata(Metadata):

    def __init__(self):
        super().__init__()
        self.compatible_platform = "agros2d"
        self.mesh_type = "triangle"

        self.nb_refinements = 1
        self.polyorder = 2
        self.adaptivity = "disabled"
        self.solver = "linear"
        self.file_suffix = '.py'


    def validate_metadata(self):
        self.validate_file_name()


class FemmMetadata(Metadata):

    def __init__(self):
        super().__init__()
        self.compatible_platform = "femm"
        self.problem_type = "electrostatic"
        self.coordinate_type = "planar"
        self.analysis_type = "steadysate"
        self.file_suffix = ".lua"

        self.frequency = 0.0
        self.unit = "m"
        self.depth = 1.0
        self.minangle=30
        self.presolven = None
        self.timestep=1e-3
        self.acsolver = 0
        self.elementsize = None
        self.smartmesh = True


    def validate_metadata(self):
        self.validate_file_name()

        if self.unit not in {"inches", "millimeters", "centimeters", "mils", "meters", "micrometers"}:
            raise ValueError(f"There is no {self.unit} unit.")