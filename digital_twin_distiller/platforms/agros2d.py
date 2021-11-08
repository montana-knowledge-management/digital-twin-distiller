import os
import subprocess
import sys
from copy import copy

from digital_twin_distiller.boundaries import BoundaryCondition, DirichletBoundaryCondition, NeumannBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import Metadata
from digital_twin_distiller.objects import Line, Node
from digital_twin_distiller.platforms.platform import Platform


class Agros2D(Platform):
    def __init__(self, m: Metadata):
        super().__init__(m)

    def __copy__(self):
        return Agros2D(copy(self.metadata))

    def comment(self, str_, nb_newline=1):
        self.file_script_handle.write(f"# {str_}")
        self.newline(nb_newline)

    def export_preamble(self):
        self.write("import agros2d as a2d")

    def export_metadata(self):
        self.write("problem = a2d.problem(clear=True)")
        self.write(f'problem.coordinate_type = "{self.metadata.coordinate_type}"')
        self.write(f'problem.mesh_type = "{self.metadata.mesh_type}"', nb_newline=2)
        if self.metadata.analysis_type == "harmonic":
            self.write(f"problem.frequency = {self.metadata.frequency}")

        self.write(f'{self.metadata.problem_type} = a2d.field("{self.metadata.problem_type}")')
        self.write(f'{self.metadata.problem_type}.analysis_type = "{self.metadata.analysis_type}"')
        self.write(f"{self.metadata.problem_type}.number_of_refinements = {self.metadata.nb_refinements}")
        self.write(f"{self.metadata.problem_type}.polynomial_order = {self.metadata.polyorder}")
        self.write(
            f'{self.metadata.problem_type}.solver = "{self.metadata.solver}"',
            nb_newline=2,
        )
        self.write("geometry = a2d.geometry", nb_newline=2)
        self.write(f'{self.metadata.problem_type}.adaptivity_type = "{self.metadata.adaptivity}"')
        if self.metadata.adaptivity != "disabled":
            self.write(
                f'{self.metadata.problem_type}.adaptivity_parameters["tolerance"] = {self.metadata.adaptivity_tol}'
            )
            self.write(
                f'{self.metadata.problem_type}.adaptivity_parameters["steps"] = {self.metadata.adaptivity_steps}'
            )

    def export_material_definition(self, mat: Material):
        mdict = {
            "magnetic_remanence_angle": mat.remanence_angle,
            "magnetic_velocity_y": mat.vy,
            "magnetic_current_density_external_real": mat.Je.real,
            "magnetic_current_density_external_imag": mat.Je.imag,
            "magnetic_permeability": mat.mu_r,
            "magnetic_conductivity": mat.conductivity,
            "magnetic_remanence": mat.remanence,
            "magnetic_velocity_angular": mat.angluar_velocity,
            "magnetic_velocity_x": mat.vx,
        }
        if self.metadata.analysis_type != "harmonic":
            mdict.pop("magnetic_current_density_external_imag")

        self.write(f'magnetic.add_material("{mat.name}", {str(mdict)})')

    def export_boundary_definition(self, boundary: BoundaryCondition):
        typename = None
        if self.metadata.problem_type == "magnetic":
            if isinstance(boundary, DirichletBoundaryCondition):
                typename = "magnetic_potential"
                A = boundary.valuedict.pop("magnetic_potential")
                boundary.valuedict["magnetic_potential_real"] = A.real
                if self.metadata.problem_type == "harmonic":
                    boundary.valuedict["magnetic_potential_imag"] = A.imag

            if isinstance(boundary, NeumannBoundaryCondition):
                typename = "magnetic_surface_current"
                A = boundary.valuedict.pop("surface_current")
                boundary.valuedict["magnetic_surface_current_real"] = A.real
                if self.metadata.problem_type == "harmonic":
                    boundary.valuedict["magnetic_surface_current_imag"] = A.imag

        self.write(
            f'{self.metadata.problem_type}.add_boundary("{boundary.name}", "{typename}", {str(boundary.valuedict)})'
        )

    def export_geometry_element(self, e, boundary=None):
        if isinstance(e, Node):
            pass

        if isinstance(e, Line):
            x0 = e.start_pt.x * self.metadata.unit
            y0 = e.start_pt.y * self.metadata.unit
            x1 = e.end_pt.x * self.metadata.unit
            y1 = e.end_pt.y * self.metadata.unit

            self.write(f"geometry.add_edge({x0}, {y0}, {x1}, {y1}", nb_newline=0)
            if boundary:
                self.write(
                    f", boundaries={{'{self.metadata.problem_type}': '{boundary}'}}",
                    nb_newline=0,
                )

            self.write(")")

    def export_block_label(self, x, y, mat: Material):
        x = self.metadata.unit * x
        y = self.metadata.unit * y
        self.write(f"geometry.add_label({x}, {y}, materials = {{'{self.metadata.problem_type}' : '{mat.name}'}})")

    def export_solving_steps(self):
        self.write("problem.solve()")
        self.write("a2d.view.zoom_best_fit()")

        self.write(f'f = open(r"{self.metadata.file_metrics_name}", "w")')

    def export_metrics(self, action, entity, variable):
        mappings = {
            "Bx": "Brx",
            "By": "Bry",
            "Br": "Brr",
            "Bz": "Brz",
            "Hx": "Hrx",
            "Hy": "Hry",
        }
        if action == "point_value":
            x = self.metadata.unit * entity[0]
            y = self.metadata.unit * entity[1]
            self.write(f'point = {self.metadata.problem_type}.local_values({x}, {y})["{mappings[variable]}"]')
            self.write(
                f'f.write("{{}}, {x}, {y}, {{}}\\n".format("{variable}", point))',
                nb_newline=2,
            )

        if action == "mesh_info":
            self.write(f"info = {self.metadata.problem_type}.solution_mesh_info()")
            self.write(f'f.write("{{}}, {{}}\\n".format("dofs", info["dofs"]))')
            self.write(f'f.write("{{}}, {{}}\\n".format("nodes", info["nodes"]))')
            self.write(f'f.write("{{}}, {{}}\\n".format("elements", info["elements"]))')

        if action == "integration":
            if self.metadata.problem_type == "magnetic":
                mapping = {"Energy": "Wm"}
                self.write(f"val={self.metadata.problem_type}.volume_integrals({entity})[{mapping[variable]!r}]")
                self.write(f'f.write("{variable}, {{}}\\n".format(val))')

    def export_closing_steps(self):
        self.write("f.close()")

    def execute(self, cleanup=False, timeout=10):
        try:
            if sys.platform == "linux":
                subprocess.run(
                    ["agros2d_solver", "-s", self.metadata.file_script_name],
                    capture_output=True,
                )
            elif sys.platform == "win32":
                subprocess.run(
                    ["Solver.exe", "-s", self.metadata.file_script_name],
                    capture_output=True,
                )

            if cleanup:
                os.remove(self.metadata.file_script_name)
                os.remove(self.metadata.file_metrics_name)
            return True
        except Exception as e:
            print(e)
            return None
