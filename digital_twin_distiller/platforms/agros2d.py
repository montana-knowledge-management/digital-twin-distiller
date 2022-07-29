import os
import subprocess
import sys
from copy import copy
from math import pi

from digital_twin_distiller.boundaries import BoundaryCondition, DirichletBoundaryCondition, NeumannBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import Metadata
from digital_twin_distiller.objects import CircleArc, Line, Node
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

        self.write("labels = []")

    def export_material_definition(self, mat: Material):
        field = self.metadata.problem_type

        if field == "electrostatic":
            mdict = {"electrostatic_charge_density": mat.Rho, "electrostatic_permittivity": mat.epsioln_r}

        if field == "magnetic":
            mdict = {
                "magnetic_remanence_angle": mat.remanence_angle,
                "magnetic_velocity_y": mat.vy,
                "magnetic_current_density_external_real": mat.Je.real,
                "magnetic_current_density_external_imag": mat.Je.imag,
                "magnetic_permeability": mat.mu_r,
                "magnetic_conductivity": mat.conductivity,
                "magnetic_remanence": mat.coercivity * 4 * pi * 1e-7,
                "magnetic_velocity_angular": mat.angluar_velocity,
                "magnetic_velocity_x": mat.vx,
            }
            if self.metadata.analysis_type != "harmonic":
                mdict.pop("magnetic_current_density_external_imag")

        if field == "heat":
            mdict = {
                "heat_density": mat.material_density,
                "heat_conductivity": mat.heat_conductivity,
                "heat_volume_heat": mat.volume_heat,
                "heat_specific_heat": mat.specific_heat,
                "heat_velocity_angular": mat.angluar_velocity,
                "heat_velocity_x": mat.vx,
                "heat_velocity_y": mat.vy,
            }

        self.write(f'{field}.add_material("{mat.name}", {str(mdict)})')

    def export_boundary_definition(self, boundary: BoundaryCondition):
        typename = None
        field = self.metadata.problem_type

        if field == "electrostatic":
            if isinstance(boundary, DirichletBoundaryCondition):
                typename = "electrostatic_potential"
                boundaryvalues = {"electrostatic_potential": boundary.valuedict["fixed_voltage"]}

            if isinstance(boundary, NeumannBoundaryCondition):
                typename = "electrostatic_surface_charge_density"
                boundaryvalues = {"electrostatic_surface_charge_density": boundary.valuedict["surface_charge_density"]}

        if field == "magnetic":
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

            boundaryvalues = boundary.valuedict.copy()

        if field == "heat":
            boundaryvalues = {
                "heat_radiation_ambient_temperature": 293.15,
                "heat_convection_external_temperature": 293.15,
                "heat_convection_heat_transfer_coefficient": 5,
                "heat_heat_flux": 0,
                "heat_radiation_emissivity": 0,
                "heat_temperature": 0,
            }
            if isinstance(boundary, DirichletBoundaryCondition):
                typename = "heat_temperature"
                boundaryvalues.clear()
                boundaryvalues["heat_temperature"] = boundary.valuedict["temperature"]

            if isinstance(boundary, NeumannBoundaryCondition):
                typename = "heat_heat_flux"  # yes, 2 x heat
                boundaryvalues.pop("heat_temperature")
                # GK: TODO: untangle the other mapping
                boundaryvalues["heat_heat_flux"] = boundary.valuedict["heat_flux"]

        self.write(f'{field}.add_boundary("{boundary.name}", "{typename}", {str(boundaryvalues)})')

    def export_geometry_element(self, e, boundary=None):
        if isinstance(e, Node):
            pass

        elif isinstance(e, Line):
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

        elif isinstance(e, CircleArc):
            x0 = e.start_pt.x * self.metadata.unit
            y0 = e.start_pt.y * self.metadata.unit
            x1 = e.end_pt.x * self.metadata.unit
            y1 = e.end_pt.y * self.metadata.unit

            xapex = e.apex_pt.x * self.metadata.unit
            yapex = e.apex_pt.y * self.metadata.unit

            self.write(f"geometry.add_edge({x0}, {y0}, {xapex}, {yapex}, angle={e.theta/2}", nb_newline=0)
            if boundary:
                self.write(f", boundaries={{'{self.metadata.problem_type}': '{boundary}'}}", nb_newline=0)
            self.write(")")

            self.write(f"geometry.add_edge({xapex}, {yapex}, {x1}, {y1}, angle={e.theta/2}", nb_newline=0)
            if boundary:
                self.write(f", boundaries={{'{self.metadata.problem_type}': '{boundary}'}}", nb_newline=0)
            self.write(")")

    def export_block_label(self, x, y, mat: Material):
        x = self.metadata.unit * x
        y = self.metadata.unit * y
        self.write(f"geometry.add_label({x}, {y}, materials = {{'{self.metadata.problem_type}' : '{mat.name}'}})")
        self.write(f"labels.append(({x}, {y}))")

    def export_solving_steps(self):
        self.write("problem.solve()")
        self.write("a2d.view.zoom_best_fit()")

        self.write(f'f = open(r"{self.metadata.file_metrics_name}", "w")')

    # TODO: check!
    def export_results(self, action, entity, variable, custom_name):
        """
        Exports the given value from the agros2d with the given coordinates.

        :param action: 'point_value', 'mesh_info', 'integration'
        :param entity:  list, looking for the closest mesh element (entity) at the given point. [x1,y1]
        :variable: W_m, exports the magnetic field energy at the given location.
        """
        mappings = {
            "Bx": "Brx",
            "By": "Bry",
            "Br": "Brr",
            "Bz": "Brz",
            "Hx": "Hrx",
            "Hy": "Hry",
            "T": "T",
            "V": "V",
            "Ex": "Ex",
            "Ey": "Ey",
        }

        custom_name_result = custom_name or variable

        field = self.metadata.problem_type
        if action == "point_value":
            x = self.metadata.unit * entity[0]
            y = self.metadata.unit * entity[1]
            self.write(f'point = {field}.local_values({x}, {y})["{mappings[variable]}"]')
            self.write(
                f'f.write("{{}}, {x}, {y}, {{}}\\n".format("{custom_name_result}", point))',
                nb_newline=2,
            )

        if action == "mesh_info":
            self.write(f"info = {self.metadata.problem_type}.solution_mesh_info()")
            self.write(f'f.write("{{}}, {{}}\\n".format("dofs", info["dofs"]))')
            self.write(f'f.write("{{}}, {{}}\\n".format("nodes", info["nodes"]))')
            self.write(f'f.write("{{}}, {{}}\\n".format("elements", info["elements"]))')

        if action == "integration":
            if field == "electrostatic":
                mapping = {"Energy": "We"}
                self.write(f"val={field}.volume_integrals({entity})[{mapping[variable]!r}]")
                self.write(f'f.write("{custom_name_result}, {{}}\\n".format(val))')

            if field == "heat":
                mapping = {"T": "T"}
                self.write(f"val={field}.volume_integrals({entity})[{mapping[variable]!r}]")
                self.write(f'f.write("{custom_name_result}, {{}}\\n".format(val))')

            if field == "magnetic":
                mapping = {"Energy": "Wm", "Torque": "Tt"}
                self.write(
                    f"selected_labels = {[(ei[0]*self.metadata.unit, ei[1]*self.metadata.unit) for ei in entity]}"
                )
                self.write(
                    f"val={field}.volume_integrals([labels.index(li) for li in selected_labels])[{mapping[variable]!r}]"
                )
                self.write(f'f.write("{custom_name_result}, {{}}\\n".format(val))')

    def export_closing_steps(self):
        self.write("f.close()")

    def execute(self, cleanup=False, timeout=10, **kwargs):
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
