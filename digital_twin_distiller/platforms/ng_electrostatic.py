from digital_twin_distiller import Material
from digital_twin_distiller.boundaries import BoundaryCondition
from digital_twin_distiller.platforms.platform import Platform


class NgElectrostatics(Platform):
    def comment(self, str_, nb_newline=1):
        pass

    def export_preamble(self):
        pass

    def export_metadata(self):
        pass

    def export_material_definition(self, mat: Material):
        pass

    def export_block_label(self, x, y, mat: Material):
        pass

    def export_boundary_definition(self, b: BoundaryCondition):
        pass

    def export_geometry_element(self, e, boundary=None):
        pass

    def export_solving_steps(self):
        pass

    def export_results(self, action, entity, variable):
        pass

    def export_closing_steps(self):
        pass

    def execute(self, cleanup=False, timeout=10):
        pass

    def __copy__(self):
        pass
