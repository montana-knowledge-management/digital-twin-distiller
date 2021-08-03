from pathlib import Path
from shutil import rmtree
from uuid import uuid4

from adze_modeler.geometry import Geometry
from adze_modeler.metadata import FemmMetadata
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot


class BaseModel:

    dir_current = Path(__file__).parent
    dir_resources = dir_current / "resources"
    dir_snapshots = dir_current / "snapshots"

    def __init__(self, X: list, i:int, j:int, exportname:str=None):
        # paths
        self.name = exportname or str(uuid4())
        self.dir_export = self.dir_snapshots / self.name
        self.file_solver_script = self.dir_export / f"P_{self.name}"
        self.file_solution = self.dir_export / f"S_{self.name}.csv"
        self.dir_snapshots.mkdir(exist_ok=True)
        self.dir_export.mkdir(exist_ok=True)

        assert len(X) == 16
        self.X = X.copy()

        # matrix indices
        self.i = i
        self.j = j

        # solver setup
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "axysimmetric"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = True
        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)
        self.g = Geometry()

        self.init_geometry()

    def define_materials(self):
        ...

    def define_boundary_conditions(self):
        ...

    def init_geometry(self):
        air = ModelPiece('air')
        air.load_piece_from_svg(self.dir_resources / "air.svg")
        air.put(495, 0)
        self.g.merge_geometry(air.geom)

        hv_coil = ModelPiece('hv_coil')
        hv_coil.load_piece_from_svg(self.dir_resources / "hv_coil.svg")
        hv_coil.put(595, 80)
        self.g.merge_geometry(hv_coil.geom)

        lv_coil = ModelPiece('lv_coil')
        lv_coil.load_piece_from_svg(self.dir_resources / "lv_coil.svg")

        offset = 0.0
        for k in range(16):
            coil = lv_coil.spawn()
            offset += 62.5 + self.X[k]
            coil.put(750, 1280-offset)
            self.g.merge_geometry(coil.geom)

        self.g.merge_lines()
        self.snapshot.add_geometry(self.g)
        self.g.export_svg(self.dir_export / 'geom.svg')

    def __call__(self, cleanup=True, timeout=1e6):
        try:
            # self.assign_blocklabels()
            # self.assign_boundary_conditions()
            # self.add_postprocesssing()
            if self.name == "dev":
                self.snapshot.export(develmode=True)
                self.snapshot.execute(cleanup=False, timeout=1e6)
            else:
                self.snapshot.export(develmode=False)
                self.snapshot.execute(cleanup=True, timeout=timeout)

            res = self.snapshot.retrive_results()

            if cleanup:
                rmtree(self.dir_export)

            if len(res) == 0:
                return None
            else:
                return res

        except Exception as e:
            print("sth went woong: ", e)
            return None




if __name__=='__main__':
    from random import seed, uniform

    seed(42)
    # X = [uniform(3, 6.25) for _ in range(16)]
    # X[0] = uniform(3, 180)
    # print(X)
    X = [5]*16
    X[0]=114
    m = BaseModel(X, i=0, j=1, exportname="dev")
    print(m())