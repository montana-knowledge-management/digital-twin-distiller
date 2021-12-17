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
        femm_metadata.depth = 1000

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        air = Material('air')
        m19 = Material('M-19')

        self.snapshot.add_material(air)
        self.snapshot.add_material(m19)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

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
