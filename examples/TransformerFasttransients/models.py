from pathlib import Path
from shutil import rmtree
from uuid import uuid4

import numpy as np

from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.geometry import Geometry
from adze_modeler.material import Material
from adze_modeler.metadata import FemmMetadata
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.objects import Rectangle


class BaseModel:

    dir_current = Path(__file__).parent
    dir_resources = dir_current / "resources"
    dir_snapshots = dir_current / "snapshots"

    def __init__(self, i: int, j: int, Ii=1.0, Ij=1.0, exportname: str = None):
        # paths
        self.name = exportname or str(uuid4())
        self.dir_export = self.dir_snapshots / self.name
        self.file_solver_script = self.dir_export / f"P_{self.name}"
        self.file_solution = self.dir_export / f"S_{self.name}.csv"
        self.dir_snapshots.mkdir(exist_ok=True)
        self.dir_export.mkdir(exist_ok=True)

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
        femm_metadata.unit = "meters"
        femm_metadata.smartmesh = True
        self.platform = Femm(femm_metadata)

        self.snapshot = Snapshot(self.platform)
        self.g = Geometry()

        self.init_geometry()
        self.define_materials()
        self.define_boundary_conditions()

    def define_materials(self):
        air = Material("air")

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)
        self.snapshot.add_boundary_condition(a0)

    def init_geometry(self):
        boundary = Rectangle(0.10, 0.3)
        boundary.put(0.075, 0, fx_point='a')



        self.g.merge_lines()
        self.snapshot.add_geometry(self.g)
        self.g.export_svg(self.dir_export / "geom.svg")

    def assign_blocklabels(self):
        ...

    def assign_boundary_conditions(self):
        ...

    def add_postprocesssing(self):
        ...

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
    Ii = 1.0
    Ij = -1.0
    m = BaseModel(i=0, j=10, Ii=Ii, Ij=Ij, exportname="dev")
    print(m())