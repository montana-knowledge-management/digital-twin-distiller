from abc import ABCMeta, abstractmethod
from pathlib import Path


class Metadata(metaclass=ABCMeta):
    def __init__(self):
        self.compatible_platform = None
        self.problem_type = None
        self.analysis_type = None
        self.coordinate_type = None
        self.mesh_type = None
        self.unit = 1.0  # 1 m
        self.precision = 1e-8
        self.frequency = 0.0

        self.file_suffix = None
        self.file_script_name = None
        self.file_metrics_name = "fem_data.csv"

    def validate_file_name(self):
        if self.file_script_name is not None:
            self.file_script_name = self.format_file_script_name() + self.file_suffix
            self.file_metrics_name = str(self.file_metrics_name)
        else:
            print("script_name is empty!")
            exit(10)

    # cut the extension of 'file_script_name' file
    def format_file_script_name(self):
        path = Path(self.file_script_name)
        return Path.joinpath(path.parent, path.stem).as_posix()

    @abstractmethod
    def validate_metadata(self):
        ...

    @abstractmethod
    def __copy__(self):
        ...


class Agros2DMetadata(Metadata):
    def __init__(self):
        super().__init__()
        self.compatible_platform = "agros2d"
        self.mesh_type = "triangle"

        self.nb_refinements = 1
        self.polyorder = 2
        self.adaptivity = "disabled"
        self.adaptivity_steps = 10
        self.adaptivity_tol = 1.0
        self.solver = "linear"
        self.file_suffix = ".py"

    def validate_metadata(self):
        self.validate_file_name()

    def __copy__(self):
        newm = Agros2DMetadata()

        newm.compatible_platform = self.compatible_platform
        newm.problem_type = self.problem_type
        newm.analysis_type = self.analysis_type
        newm.coordinate_type = self.coordinate_type
        newm.mesh_type = self.mesh_type
        newm.unit = self.unit
        newm.precision = self.precision
        newm.file_suffix = self.file_suffix
        newm.file_script_name = self.file_script_name
        newm.file_metrics_name = self.file_metrics_name

        newm.compatible_platform = self.compatible_platform
        newm.mesh_type = self.mesh_type
        newm.nb_refinements = self.nb_refinements
        newm.polyorder = self.polyorder
        newm.adaptivity = self.adaptivity
        newm.adaptivity_tol = self.adaptivity_tol
        newm.adaptivity_steps = self.adaptivity_steps
        newm.solver = self.solver
        newm.file_suffix = self.file_suffix
        return newm


class FemmMetadata(Metadata):
    def __init__(self):
        super().__init__()
        self.compatible_platform = "femm"
        self.problem_type = "electrostatic"
        self.coordinate_type = "planar"
        self.analysis_type = "steadysate"
        self.file_suffix = ".lua"

        self.unit = "meters"
        self.depth = 1.0
        self.minangle = 30
        self.presolven = None
        self.timestep = 1e-3
        self.acsolver = 0
        self.elementsize = None
        self.smartmesh = True

    def validate_metadata(self):
        self.validate_file_name()

        if self.unit not in {
            "inches",
            "millimeters",
            "centimeters",
            "mils",
            "meters",
            "micrometers",
        }:
            raise ValueError(f"There is no {self.unit} unit.")

    def __copy__(self):
        newm = FemmMetadata()

        newm.compatible_platform = self.compatible_platform
        newm.problem_type = self.problem_type
        newm.analysis_type = self.analysis_type
        newm.coordinate_type = self.coordinate_type
        newm.mesh_type = self.mesh_type
        newm.unit = self.unit
        newm.precision = self.precision
        newm.file_suffix = self.file_suffix
        newm.file_script_name = self.file_script_name
        newm.file_metrics_name = self.file_metrics_name

        newm.compatible_platform = self.compatible_platform
        newm.problem_type = self.problem_type
        newm.coordinate_type = self.coordinate_type
        newm.analysis_type = self.analysis_type
        newm.file_suffix = self.file_suffix
        newm.frequency = self.frequency
        newm.unit = self.unit
        newm.depth = self.depth
        newm.minangle = self.minangle
        newm.presolven = self.presolven
        newm.timestep = self.timestep
        newm.acsolver = self.acsolver
        newm.elementsize = self.elementsize
        newm.smartmesh = self.smartmesh
        return newm


class NgSolveMetadata(Metadata):
    def __init__(self):
        super().__init__()
        self.compatible_platform = "ngsolve"
        self.coordinate_type = "planar"
        self.problem_type = "electrostatic"
        self.file_suffix = ".py"
        self.depth = 1

    def validate_metadata(self):
        self.validate_file_name()

    def __copy__(self):
        ...


class NgElectrostaticMetadata(NgSolveMetadata):
    def __init__(self):
        super().__init__()
        self.compatible_platform = "ngsolve"
        self.coordinate_type = "planar"
        self.problem_type = "electrostatic"
        self.file_suffix = ".py"
        self.depth = 1

    def validate_metadata(self):
        self.validate_file_name()

    def __copy__(self):
        ...
