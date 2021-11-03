from copy import copy
from pathlib import Path
from shutil import rmtree
from uuid import uuid4

from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.geometry import Geometry
from adze_modeler.material import Material
from adze_modeler.metadata import Agros2DMetadata, FemmMetadata
from adze_modeler.objects import Rectangle
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot


class BaseModel:

    dir_current = Path(__file__).parent
    dir_resources = dir_current / "resources"
    dir_snapshots = dir_current / "snapshots"

    def __init__(
        self, i: int = -1, j: int = -1, Ii=1.0, Ij=1.0, exportname: str = None
    ):
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
        if i == j:
            self.Ij = 0.0
            self.j = -1

        # solver setup
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "meters"
        femm_metadata.smartmesh = True
        femm_metadata.frequency = 10e6
        femm_metadata.precision = 1e-12
        self.platform = Femm(femm_metadata)

        # agros_metadata = Agros2DMetadata()
        # agros_metadata.file_script_name = self.file_solver_script
        # agros_metadata.file_metrics_name = self.file_solution
        # agros_metadata.problem_type = "magnetic"
        # agros_metadata.coordinate_type = "planar"
        # agros_metadata.analysis_type = "steadystate"
        # agros_metadata.unit = 1
        # agros_metadata.nb_refinements = 0
        # agros_metadata.polyorder = 2
        # agros_metadata.adaptivity = "hp-adaptivity"
        # agros_metadata.adaptivity_tol = 0.001
        # agros_metadata.adaptivity_steps = 30
        # self.platform = Agros2D(agros_metadata)

        self.snapshot = Snapshot(self.platform)
        self.g = Geometry()

        # reusable geometry elements
        self.r_window = Rectangle(width=0.1, height=0.3)
        self.r_window.put(0.075, 0.075, fx_point="a")

        self.coil_i = Rectangle(width=0.004, height=0.004)
        self.coil_j = Rectangle(width=0.004, height=0.004)

        self.init_geometry()
        self.define_materials()
        self.define_boundary_conditions()

    def define_materials(self):
        air = Material("air")

        coil_area = self.coil_i.width * self.coil_i.height
        Ii = Material("Ii")
        Ii.Je = self.Ii / coil_area

        Ij = Material("Ij")
        Ij.Je = self.Ij / coil_area

        self.snapshot.add_material(air)
        self.snapshot.add_material(Ii)
        if self.j != -1 and self.i != self.j:
            self.snapshot.add_material(Ij)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition(
            "a0", field_type="magnetic", magnetic_potential=0.0
        )
        self.snapshot.add_boundary_condition(a0)

    def init_geometry(self):
        self.g.add_rectangle(self.r_window)

        self.coil_i.put(0.075 + 0.018, 0.075 + 0.3 - 0.03 - 0.008 * self.i)
        self.g.add_rectangle(self.coil_i)

        if self.j != -1 and self.i != self.j:
            self.coil_j.put(0.075 + 0.018, 0.075 + 0.3 - 0.03 - 0.008 * self.j)
            self.g.add_rectangle(self.coil_j)

        self.g.merge_lines()
        self.snapshot.add_geometry(self.g)
        self.g.export_svg(self.dir_export / "geom.svg")

    def assign_blocklabels(self):
        self.snapshot.assign_material(*self.r_window.cp, "air")
        self.snapshot.assign_material(*self.coil_i.cp, "Ii")
        if self.j != -1 and self.i != self.j:
            self.snapshot.assign_material(*self.coil_j.cp, "Ij")

    def assign_boundary_conditions(self):
        left = (self.r_window.a + self.r_window.d) * 0.5
        right = (self.r_window.b + self.r_window.c) * 0.5
        upper = (self.r_window.d + self.r_window.c) * 0.5
        lower = (self.r_window.a + self.r_window.b) * 0.5

        self.snapshot.assign_boundary_condition(*left, "a0")
        self.snapshot.assign_boundary_condition(*right, "a0")
        self.snapshot.assign_boundary_condition(*upper, "a0")
        self.snapshot.assign_boundary_condition(*lower, "a0")

    def add_postprocesssing(self):
        window = list(self.r_window.cp)

        self.snapshot.add_postprocessing("mesh_info", None, None)

        # use this line if platform==agros2d:
        self.snapshot.add_postprocessing("integration", [window], "Energy")

        # use this line if platform==agros2d
        # self.snapshot.add_postprocessing("integration", [0], "Energy")

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


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np

    Ii = 1.0
    Ij = 1.0

    # simple model for 1 calculation
    # m = BaseModel(i=0, j=5, Ii=Ii, Ij=Ij, exportname="dummyname")
    # E = m(cleanup=False).get('Energy')
    # print(2*2*E/(Ii**2))

    # Calculation of the energy matrix
    # E = np.zeros((30, 30))
    # for i in range(30):
    #     for j in range(i, 30):
    #         m = BaseModel(i=i, j=j, Ii=Ii, Ij=Ij)
    #         E[i, j] = 2 * m(cleanup=True).get('Energy')
    #         E[j, i] = E[i, j]
    #         print(f'{E[i,j]:.2e}', end=' ')
    #     print()

    # np.savetxt(BaseModel.dir_resources / "Em.txt", E)
    # plt.matshow(E)
    # # plt.savefig(BaseModel.dir_current / "docs/media/Ematrix.png")
    # plt.show()

    # Calculation of the inductances based on the energy matrix
    E = np.loadtxt(BaseModel.dir_resources / "Em.txt")
    M = np.zeros_like(E)

    # calculate the self inductanes first
    for i in range(30):
        Lii = 2 * E[i, i] / (Ii ** 2)
        M[i, i] = Lii

    # calculate the mutual inductances
    for j in range(30):
        for k in range(j + 1, 30):
            M[j, k] = (E[j, k] - (E[j, j] + E[k, k])) / (Ii * Ij)
            # M[j, k] = 1/(Ii*Ij)*(E[j, k] - 0.5*Ii*M[j, j]-0.5*Ij*M[k, k])
            M[k, j] = M[j, k]

    np.savetxt(BaseModel.dir_resources / "M.txt", M)

    plt.matshow(M)
    plt.colorbar()
    plt.savefig(
        BaseModel.dir_current / "docs/media/inductance_matrix.png",
        bbox_inches="tight",
        dpi=550,
    )
    plt.show()

    fix, axes = plt.subplots(3, sharex=True, figsize=(10, 7))
    coils = range(1, 31)

    axes[0].bar(coils, M[0, :], color="b", label="Coil 1")
    axes[0].grid()
    axes[0].legend()
    axes[0].set_ylabel("Inductance")

    axes[1].bar(coils, M[14, :], color="r", label="Coil 15")
    axes[1].grid()
    axes[1].legend()
    axes[1].set_ylabel("Inductance")

    axes[2].bar(coils, M[29, :], color="g", label="Coil 30")
    axes[2].grid()
    axes[2].legend()
    axes[2].set_ylabel("Inductance")

    plt.xlabel("Coils")
    plt.xticks(coils)
    plt.savefig(
        BaseModel.dir_current / "docs/media/inductance_change.png",
        bbox_inches="tight",
        dpi=550,
    )
    plt.show()

    # Mesh density vs self inductance
    # nb_nodes = []
    # energies = []
    # mesh_sizes = np.linspace(0.5e-3, 5e-3, 21)
    # for mesh_size_i in mesh_sizes:
    #     m = BaseModel(i=0, Ii=Ii, Ij=Ij, exportname="mesh")
    #     m.snapshot.materials['air'].meshsize = mesh_size_i
    #     res = m(cleanup=False)
    #     nb_nodes.append(res['nodes'])
    #     energies.append(res['Energy'])
    #     print(mesh_size_i, nb_nodes[-1], res['Energy'])

    # plt.figure(figsize=(7,4))
    # plt.plot(nb_nodes, energies, 'b-o', lw=2)
    # plt.xlabel('Number of nodes')
    # plt.ylabel(r'Self inductance (L$_{11}$)')
    # plt.grid()
    # plt.savefig(BaseModel.dir_current / "docs/media/meshsize_vs_selfinductance.png")
    # plt.show()
