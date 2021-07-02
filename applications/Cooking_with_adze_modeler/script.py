from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.material import Material
from adze_modeler.model import Model
from adze_modeler.geometry import Geometry
from random import uniform
from pathlib import Path
from adze_modeler.metadata import FemmMetadata, Agros2DMetadata
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from numpy import linspace, meshgrid
from uuid import uuid4

from statistics import stdev

class CoilProblem(Model):

    def __init__(self):
        super().__init__()
        self.export_location = Path(__file__).parent / "snapshots"

    def build(self, p: dict, cleanup=False):
        """
        Let's write the recipie
        """

        self.model_id = str(uuid4())

        # 1. Create a snapshot object
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "axisymmetric"
        femm_metadata.file_script_name = self.export_location / self.model_id / "femm_solver_script"
        femm_metadata.file_metrics_name = self.export_location / self.model_id / "femm_solution.csv"
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False

        agros_metadata = Agros2DMetadata()
        agros_metadata.file_script_name = self.export_location / self.model_id / "agros_solver_script"
        agros_metadata.file_metrics_name = self.export_location / self.model_id /  "agros_solution.csv"
        agros_metadata.problem_type = "magnetic"
        agros_metadata.coordinate_type = "axisymmetric"
        agros_metadata.analysis_type = "steadystate"
        agros_metadata.unit = 1e-3
        # agros_metadata.nb_refinements = 2
        # agros_metadata.adaptivity = "hp-adaptivity"

        platform_femm = Femm(femm_metadata)
        platform_agros = Agros2D(agros_metadata)
        self.snapshot = Snapshot(platform_femm)

        # 1.1 Define materials
        exctitation = Material("J+")
        exctitation.Je = 2e6
        air = Material("air")
        core = Material("core")
        core.mu_r = 100
        core.meshsize = 0.1
        self.snapshot.add_material(exctitation)
        self.snapshot.add_material(air)
        self.snapshot.add_material(core)

        self.snapshot.assign_material(3, 0, name="core")
        self.snapshot.assign_material(30, 30, name="air")

        # 1.2 Define Boundary Conditions
        b1 = DirichletBoundaryCondition(name="a0", field_type=self.snapshot.platform.metadata.problem_type,
                                        magnetic_potential=0.0)
        # print(b1)



        # 2. Create a geometry using the modelpieces from the ingredients
        if isinstance(p, dict):
            N = len(p.keys())
        else:
            N = len(p)
        h = 1.5
        geom = Geometry()
        geom.merge_geometry(self.ingredients.get("bound").geom)
        geom.merge_geometry(self.ingredients.get("core").geom)
        # for i, ri in enumerate(p.values()):
        for i, ri in enumerate(p):
            coil = self.ingredients.get('coil').spawn()
            offsetz = i * h - N * h / 2
            coil.put(ri, offsetz)
            geom.merge_geometry(coil.geom)
            self.snapshot.assign_material(ri + h/2, offsetz + h/2, name="J+")

        geom.generate_intersections()
        self.snapshot.add_geometry(geom)

        self.snapshot.add_boundary_condition(b1)
        self.snapshot.assign_boundary_condition(0, 0, "a0")
        self.snapshot.assign_boundary_condition(0, -20, "a0")
        self.snapshot.assign_boundary_condition(0, 20, "a0")
        self.snapshot.assign_boundary_condition(35, 35, "a0")
        self.snapshot.assign_boundary_condition(35, -35, "a0")
        self.snapshot.assign_boundary_condition(70, 0, "a0")

        self.model_path = self.export_location / f"{self.model_id}"
        geom.export_geom(self.model_path / f"{self.model_id}.svg")

        # aadding postprocessing steps
        # adding measurement points to the core
        Nx = 10
        Ny = 10
        px = linspace(0, 5, Nx)
        py = linspace(-5, 5, Ny)
        xv, yv = meshgrid(px, py, sparse=False, indexing='xy')

        for i in range(Nx):
            for j in range(Ny):
                eval_point = (xv[j, i], yv[j, i])
                self.snapshot.add_postprocessing('point_value', eval_point, 'Bz')

        self.snapshot.add_postprocessing("mesh_info", None, None)
        self.snapshot.export()
        self.snapshot.execute()
        res = self.snapshot.retrive_results()
        if cleanup:
            self.cleanup()

        B0 = 2e-3
        F1 = max(map(lambda pointvalue: abs(pointvalue[2] - B0), res['Bz']))
        # F1 = stdev([pointvalue[2] for pointvalue in res['Bz']])

        # print(self.snapshot.platform.metadata.compatible_platform)
        # print("number of nodes:", res['nodes'])
        # print("number of elements:", res['elements'])
        # print("F1:", F1)
        # print('-'*50 + '\n')
        print(F1)
        return F1



# Ingredients
mp_bound = ModelPiece("bound")
mp_bound.load_piece_from_svg('problem_boundary.svg')

mp_core = ModelPiece("core")
mp_core.load_piece_from_svg('core.svg')
mp_core.put(0, -5)

mp_coil = ModelPiece("coil")
mp_coil.load_piece_from_svg('coil.svg')


material_air = Material("air")
material_core = Material("core")
material_core.mu_r = 1200
material_J = Material("J+")
material_J.Je = 2e6


m = CoilProblem()
m.add_piece(mp_bound)
m.add_piece(mp_coil)
m.add_piece(mp_core)

N =20
r0 = [uniform(5, 50) for _ in range(N)]
parameters = {f'r{i}': ri for i, ri in enumerate(r0)}
for namei, ri in parameters.items():
    print(namei, ri)

print(m.build(r0, cleanup=False))

M = 20
for i in range(M):
    r0 = [uniform(5, 50) for _ in range(N)]
    m.build(r0, cleanup=False)

# from scipy.optimize import minimize
# r0 = [uniform(5, 50) for _ in range(N)]
# res = minimize(m.build, r0, method="SLSQP", bounds=[(5, 50) for _ in range(N)], tol=0.0005, options={'maxiter':50})
# print(res)
#
# print(m.build(res.x, cleanup=False))