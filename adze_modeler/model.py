from abc import abstractmethod
from adze_modeler import objects as obj
from adze_modeler.geometry import Geometry
from adze_modeler.snapshot import Snapshot
from shutil import rmtree
from uuid import uuid4
from pathlib import Path
import sys

class BaseModel:

    def __init__(self, exportname=None): 
        self.name = exportname or str(uuid4())

        self.dir_current = Path(sys.modules[self.__module__].__file__).parent
        self.dir_resources = self.dir_current / "resources"
        self.dir_snapshots = self.dir_current / "snapshots"
        self.dir_media = self.dir_current / "docs" / "media"
        self.dir_data = self.dir_current / "data"
        self.dir_export = self.dir_snapshots / self.name


        self.file_solver_script = self.dir_export / f"P_{self.name}"
        self.file_solution = self.dir_export / f"S_{self.name}.csv"

        self.snapshot: Snapshot = None
        self.geom = Geometry()

        self._init_directories()

    def add_line(self, x0, y0, x1, y1):
        self.g.add_line(obj.Line(obj.Node(x0, y0), obj.Node(x1, y1)))

    def _init_directories(self):
        self.dir_resources.mkdir(exist_ok=True)
        self.dir_snapshots.mkdir(exist_ok=True)
        self.dir_media.mkdir(exist_ok=True, parents=True)
        self.dir_data.mkdir(exist_ok=True)
        self.dir_snapshots.mkdir(exist_ok=True)
        self.dir_export.mkdir(exist_ok=True)

    @abstractmethod
    def setup_solver(self):
        ...

    @abstractmethod
    def add_postprocessing(self):
        ...

    @abstractmethod
    def define_materials(self):
        ...

    @abstractmethod
    def assign_materials(self):
        ...
        
    @abstractmethod
    def define_boundary_conditions(self):
        ...

    @abstractmethod
    def assign_boundary_conditions(self):
        ...

    @abstractmethod
    def build_geometry(self):
       ... 

    def __call__(self, cleanup=True, devmode=False, timeout=100):
        try:
            self.setup_solver()
            self.define_materials()
            self.define_boundary_conditions()
            self.build_geometry()
            self.assign_materials()
            self.assign_boundary_conditions()
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
            return None
