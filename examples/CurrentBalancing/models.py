from pathlib import Path
from shutil import rmtree
from uuid import uuid4

import numpy as np

from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.geometry import Geometry
from adze_modeler.material import Material
from adze_modeler.metadata import FemmMetadata, Agros2DMetadata
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot


class BaseModel:

    dir_current = Path(__file__).parent
    dir_resources = dir_current / "resources"
    dir_snapshots = dir_current / "snapshots"

    def __init__(self, X: list, i: int, j: int, Ii=1.0, Ij=1.0, exportname: str = None):
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
        self.Ii = Ii
        self.j = j
        self.Ij = Ij

        # solver setup
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "axisymmetric"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = True
        self.platform = Femm(femm_metadata)

        self.snapshot = Snapshot(self.platform)
        self.g = Geometry()

        self.init_geometry()
        self.define_materials()
        self.define_boundary_conditions()

    def define_materials(self):
        air = Material("air")
        hv = Material("Primary")
        hv.Je = 30e3 / (30 * 1100) * 1e6
        hv.Je = 5 / (30 * 1100) * 1e6

        ii = Material("Ii")
        ii.Je = self.Ii / (33 * 62.5) * 1e6

        ij = Material("Ij")
        ij.Je = self.Ij / (33 * 62.5) * 1e6

        self.snapshot.add_material(air)
        self.snapshot.add_material(hv)
        self.snapshot.add_material(ii)
        self.snapshot.add_material(ij)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)
        self.snapshot.add_boundary_condition(a0)

    def init_geometry(self):
        air = ModelPiece("air")
        air.load_piece_from_svg(self.dir_resources / "air.svg")
        air.put(495, 0)
        self.g.merge_geometry(air.geom)

        hv_coil = ModelPiece("hv_coil")
        hv_coil.load_piece_from_svg(self.dir_resources / "hv_coil.svg")
        hv_coil.put(595, 80)
        self.g.merge_geometry(hv_coil.geom)

        lv_coil = ModelPiece("lv_coil")
        lv_coil.load_piece_from_svg(self.dir_resources / "lv_coil.svg")

        offset = 0.0
        for k in range(16):
            coil = lv_coil.spawn()
            offset += 62.5 + self.X[k]
            coil.put(750, 1280 - offset)
            self.g.merge_geometry(coil.geom)

        self.g.merge_lines()
        self.snapshot.add_geometry(self.g)
        self.g.export_svg(self.dir_export / "geom.svg")

    def assign_blocklabels(self):
        self.snapshot.assign_material(605, 699, "Primary")
        self.snapshot.assign_material(523, 1036, "air")
        offset = 0.0
        for k in range(16):
            offset += 62.5 + self.X[k]
            if k == self.i:
                self.snapshot.assign_material(750 + 33.0 / 2, 1280 - offset + 62.5 / 2, "Ii")
            elif k == self.j:
                self.snapshot.assign_material(750 + 33.0 / 2, 1280 - offset + 62.5 / 2, "Ij")
            else:
                self.snapshot.assign_material(750 + 33.0 / 2, 1280 - offset + 62.5 / 2, "air")

    def assign_boundary_conditions(self):
        self.snapshot.assign_boundary_condition(493, 695, "a0")
        self.snapshot.assign_boundary_condition(758, 1280, "a0")
        self.snapshot.assign_boundary_condition(758, 0, "a0")
        self.snapshot.assign_boundary_condition(996, 644, "a0")

    def add_postprocesssing(self):
        entities = [(605, 699), (523, 1036)]
        offset = 0.0
        for k in range(16):
            offset += 62.5 + self.X[k]
            entities.append((750 + 33.0 / 2, 1280 - offset + 62.5 / 2))
        self.snapshot.add_postprocessing("integration", entities, "Energy")

    def __call__(self, cleanup=True, timeout=1e6):
        try:
            self.assign_blocklabels()
            self.assign_boundary_conditions()
            self.add_postprocesssing()
            if self.name == "dev":
                self.snapshot.export(develmode=True)
                self.snapshot.execute(cleanup=False, timeout=1e6)
            else:
                self.snapshot.export(develmode=False)
                self.snapshot.execute(cleanup=cleanup, timeout=timeout)

            res = self.snapshot.retrive_results()

            if cleanup:
                rmtree(self.dir_export)

            if len(res) == 0:
                return None
            else:
                return res.get("Energy")

        except Exception as e:
            print("sth went woong: ", e)
            return None


if __name__ == "__main__":
    from random import seed, uniform
    from numpy import zeros
    import matplotlib.pyplot as plt

    seed(42)
    # X = [uniform(3, 6.25) for _ in range(16)]
    # X[0] = uniform(3, 180)
    # print(X)
    X = [5] * 16
    X[0] = 114

    E = zeros((16, 16))
    M = zeros((16, 16))

    # m = BaseModel(X, i=1, j=10, exportname="dev")
    for i in range(0, 16):
        for j in range(i, 16):
            m = BaseModel(X, i=i, j=j)
            Eij = m(cleanup=True)
            E[i, j] = Eij
            E[j, i] = Eij
            print(f"{Eij:.3e}", end=" ")
        print()

    print(E)
    np.savetxt(BaseModel.dir_resources / "Em.txt", E)
    plt.matshow(E)
    plt.savefig(BaseModel.dir_current / "docs/media/Ematrix.png")
    plt.show()
