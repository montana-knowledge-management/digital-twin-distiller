from digital_twin_distiller.boundaries import DirichletBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot
from digital_twin_distiller import ModelPiece
import numpy as np

ModelDir.set_base(__file__)

class FrozenPermeability(BaseModel):
    """docstring for FrozenPermeability"""
    def __init__(self, **kwargs):
        super(FrozenPermeability, self).__init__(**kwargs)
        self._init_directories()

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = True
        femm_metadata.depth = 1000 # TODO: Correct this value

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        # creating default materials
        air = Material('air')
        steel_m19 = Material('M-19')
        steel_1018 = Material('1018')
        magnet = Material('Magnet')
        coil = Material('coil')

        # default material configuration
        steel_m19.thickness = 0.635
        steel_m19.lamination_type = "inplane"
        steel_m19.fill_factor = 0.98
        steel_m19.conductivity = 1.9e6
        steel_m19.b = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
                       0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55,
                       1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30]
        steel_m19.h = [0.0, 15.120714, 22.718292, 27.842733, 31.871434, 35.365044, 38.600588, 41.736202, 44.873979,
                       48.087807, 51.437236, 54.975221, 58.752993, 62.823644, 67.245285, 72.084406, 77.420100,
                       83.350021, 89.999612, 97.537353, 106.201406, 116.348464, 128.547329, 143.765431, 163.754169,
                       191.868158, 234.833507, 306.509769, 435.255202, 674.911968, 1108.325569, 1813.085468,
                       2801.217421, 4053.653117, 5591.106890, 7448.318413, 9708.815670, 12486.931615, 16041.483644,
                       21249.420624, 31313.495878, 53589.446877, 88477.484601, 124329.410540, 159968.569300,
                       197751.604272, 234024.751347]

        steel_1018.conductivity = 5.8e6
        steel_1018.b = [0.0, 0.2503, 0.925, 1.250, 1.390, 1.525, 1.71, 1.87, 1.955, 2.02, 2.11, 2.225, 2.43]
        steel_1018.h = [0.0, 238.7325, 795.775, 1591.55, 2387.325, 3978.875, 7957.75, 15915.5, 23873.25, 39788.75,
                        79577.5, 159155, 318310]
        steel_1018.phi_hmax = 20

        magnet.mu_r = 1.04496
        magnet.coercivity = 891000

        coil.diameter = 0.643808017739015
        coil.lamination_type = "magnetwire"
        coil.conductivity = 58e6



        self.snapshot.add_material(air)
        self.snapshot.add_material(steel_m19)
        self.snapshot.add_material(steel_1018)
        self.snapshot.add_material(magnet)


    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

        self.assign_boundary_arc(0,     53, "a0")
        self.assign_boundary_arc(0,    -53, "a0")
        self.assign_boundary_arc(0,  12.35, "a0")
        self.assign_boundary_arc(0, -12.35, "a0")

    def add_postprocessing(self):
        r = 4 # amend
        phi = np.linspace(0, 2*np.pi, 101)
        for phi_i in phi:
            self.snapshot.add_postprocessing("point_value", (r*np.cos(phi_i), r*np.sin(phi_i)), "Bx")
            self.snapshot.add_postprocessing("point_value", (r*np.cos(phi_i), r*np.sin(phi_i)), "By")

    def build_geometry(self):
        stator = ModelPiece('stator')
        stator.load_piece_from_dxf(ModelDir.RESOURCES / 'stator.dxf')

        rotor = ModelPiece('rotor')
        rotor.load_piece_from_dxf(ModelDir.RESOURCES / 'rotor.dxf')
        rotor.rotate(alpha=0.0)

        self.geom.merge_geometry(stator.geom)
        self.geom.merge_geometry(rotor.geom)

        self.snapshot.add_geometry(self.geom)


if __name__ == "__main__":
    m = FrozenPermeability(exportname="dev")
    print(m(cleanup=False, devmode=True))
