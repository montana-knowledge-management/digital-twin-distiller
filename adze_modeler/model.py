from abc import ABCMeta, abstractmethod
from adze_modeler import objects as obj
from adze_modeler.geometry import Geometry
from adze_modeler.snapshot import Snapshot
from shutil import rmtree
from uuid import uuid4
from pathlib import Path
import sys
import traceback


class BaseModel(metaclass=ABCMeta):
    """
    This abstract class servers as a baseline to describe an adze-model. It also provides automatic
    path creation, model uilding, execution, results extraction and cleanup.


    """

    def __init__(self, exportname: str = None):
        """
        This function sets the paths and file names.

        Parameters:
            exportname: A specific name for a Model instance instead of a random generated string.

        """
        self.name = exportname or str(uuid4())

        self.dir_current = Path(sys.modules[self.__module__].__file__).parent
        self.dir_resources = self.dir_current / "resources"
        self.dir_snapshots = self.dir_current / "snapshots"
        self.dir_media = self.dir_current / "media"
        self.dir_data = self.dir_current / "data"
        self.dir_export = self.dir_snapshots / self.name

        self.file_solver_script = self.dir_export / f"P_{self.name}"
        self.file_solution = self.dir_export / f"S_{self.name}.csv"

        self.snapshot: Snapshot = None
        self.geom = Geometry()

        self.label_queue = []
        self.boundary_queue = []
        self.boundary_arc_queue = []

    def add_line(self, x0: float, y0: float, x1: float, y1: float):
        """
        Conviniently add a line to the `geom` attribute.

        Parameters:
            x0: x coordinate of the starting point
            y0: y coordinate of the starting point
            x1: x coordinate of the end point
            y1: y coordinate of the end point
        """
        self.geom.add_line(obj.Line(obj.Node(x0, y0), obj.Node(x1, y1)))

    def _init_directories(self):
        """
        Make the specified directories. This function will not raise an exception when a directory is already present.
        """
        self.dir_resources.mkdir(exist_ok=True)
        self.dir_snapshots.mkdir(exist_ok=True)
        self.dir_media.mkdir(exist_ok=True, parents=True)
        self.dir_data.mkdir(exist_ok=True)
        self.dir_snapshots.mkdir(exist_ok=True)
        self.dir_export.mkdir(exist_ok=True)

        self.file_solver_script = self.file_solver_script.absolute()
        self.file_solution = self.file_solution.absolute()

    def _assign_materials(self):
        """
        Iterate over the `label_queue` and assign the materials in the `snapshot` attribute.
        """

        while self.label_queue:
            self.snapshot.assign_material(*self.label_queue.pop())

    def _assign_boundary_conditions(self):
        """
        Iterate over the boundary conditions and assign them in the `snapshot` attribute.
        """

        while self.boundary_queue:
            self.snapshot.assign_boundary_condition(*self.boundary_queue.pop())

        while self.boundary_arc_queue:
            self.snapshot.assign_arc_boundary_condition(*self.boundary_arc_queue.pop())

    @abstractmethod
    def setup_solver(self):
        """
        In this function you have to initialize the `snapshot` variable with the FEM solver of your choice. You have to 
        create a metadata object that holds the setups for the particular problem. Then using this object you have to initialize a
        platform object that fits the metadata (Agros2DMetadata for Agros2D platform). Then use this platform object to initialize
        the `snapshot` variable.
        """
        ...

    @abstractmethod
    def add_postprocessing(self):
        """
        Use the `snapshot.add_postprocessing(action, entity, variable)` method to append postprocessing steps.
        """
        ...

    @abstractmethod
    def define_materials(self):
        """
        Define and add the model specific materials to the `snapshot` variable.
        """
        ...

    @abstractmethod
    def define_boundary_conditions(self):
        """
        Define and add the model specific boundary conditions to the `snapshot` variable.
        """
        ...

    @abstractmethod
    def build_geometry(self):
        """
        This is function is responsible for building the geometry. After the building is done, use the `snapshot.add_geometry` method to merge 
        your geometry into the snapshot.
        """
        ...

    def __call__(self, cleanup=True, devmode=False, timeout=100):
        """
        Calling a Model instance will execute the solution steps.

        Parameters:
            cleanup: If `True`, the the dir_export and all of its content will be deleted after the execution.
            devmode: If `True`, the solver will open the script but won't execute it. If `True`, then `cleanup`
                     is automatically set to `False`.
        """
        try:
            self.setup_solver()
            self.define_materials()
            self.define_boundary_conditions()
            self.build_geometry()
            self._assign_materials()
            self._assign_boundary_conditions()
            self.add_postprocessing()

            if devmode:
                self.snapshot.export(develmode=True)
                self.snapshot.execute(cleanup=False, timeout=1e7)
            else:
                self.snapshot.export(develmode=False)
                self.snapshot.execute(cleanup=False, timeout=timeout)

            res = self.snapshot.retrive_results()

            if cleanup:
                rmtree(self.dir_export)

            if len(res) == 0:
                return None
            else:
                return res

        except Exception as e:
            print("something went wrong: ", e)
            print(traceback.format_exc())
            return None
