from shutil import rmtree
from adze_modeler.snapshot import Snapshot
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.model import BaseModel
from adze_modeler.metadata import FemmMetadata
from adze_modeler.platforms.femm import Femm
from adze_modeler.material import Material
from adze_modeler.boundaries import AntiPeriodicBoundaryCondition, DirichletBoundaryCondition

class PriusMotor(BaseModel):
    """docstring for PriusMotor"""
    def __init__(self, exportname=None):
        super(PriusMotor, self).__init__(exportname=exportname)

    
    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = True
        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)
    
    def add_postprocessing(self):
        ...
    
    def define_materials(self):
        air = Material('air')
        core = Material('core')
        magnet = Material('magnet')

        self.snapshot.add_material(air)
        self.snapshot.add_material(core)
        self.snapshot.add_material(magnet)
        
    
    def assign_materials(self):
        ...
        
    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type='magnetic', magnetic_potential=0.0)
        pb1 = AntiPeriodicBoundaryCondition('PB1', field_type='magnetic')
        pb2 = AntiPeriodicBoundaryCondition('PB2', field_type='magnetic')
        pb3 = AntiPeriodicBoundaryCondition('PB3', field_type='magnetic')

        self.snapshot.add_boundary_condition(a0)
        self.snapshot.add_boundary_condition(pb1)
        self.snapshot.add_boundary_condition(pb2)
        self.snapshot.add_boundary_condition(pb3)

    def assign_boundary_conditions(self):
        ...

    def build_geometry(self):
        rotor_slice = ModelPiece('rotor_slice')
        rotor_slice.load_piece_from_dxf(self.dir_resources / "rotor_slice.dxf")
        rotor_slice.put(-30.6912, 51.27)

        stator_slice = ModelPiece('rotor_slice')
        stator_slice.load_piece_from_dxf(self.dir_resources / "stator_slice.dxf")
        stator_slice.put(-51.5168, 75)


        self.geom.merge_geometry(rotor_slice.geom)
        self.geom.merge_geometry(stator_slice.geom)
        self.snapshot.add_geometry(self.geom)



if __name__ == "__main__":
    m = PriusMotor()
    print(m(cleanup=False, devmode=True))