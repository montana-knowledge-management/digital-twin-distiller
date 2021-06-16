"""
to execture a script use:
    agros2d_solver -s script.py

execute a script with gui:
    agros2d -s script.py
"""


class Agros2DWrapper:
    def __init__(self):
        self.fields = {
            "electrostatic": {},
            "magnetic": {},
            "current": {},
            "heat": {},
        }
        self.coordinate_type = "planar"
        self.mesh_type = "triangle"
        self.boundary_conditions = {}

    def add_field_electrostatic(self, **kwargs):
        """
        :param analysis: 'steady'
        :param solver: 'linear'
        :param matrix_solver: 'umfpack', 'mumps', 'paraluation_it', 'paraluation_amg', 'mums_out'
        :param nb_refinements: Greater than 0 integer.
        :param nb_polyorder: Greater than 0 integer.
        :param adaptivity: 'disabled', 'h', 'p', 'hp'
        """

        analysis = kwargs.get("analysis", "steady")
        solver = kwargs.get("solver", "linear")
        matrix_solver = kwargs.get("matrix_solver", "umfpack")
        nb_refinements = kwargs.get("nb_refinements", 1)
        nb_polyorder = kwargs.get("nb_polyorder", 2)
        adaptivity = kwargs.get("adaptivity", "disabled")

        if analysis != "steady":
            raise ValueError(f"Electric field can can be: 'steady'. Got '{analysis}' instead.")

        if solver != "linear":
            raise ValueError(f"Only 'linear' solver is supported.")

        if matrix_solver not in {"umfpack", "mumps", "paraluation_it", "paraluation_amg", "mums_out"}:
            raise ValueError(f'Imporper matrix_solver type: "{matrix_solver}"')

        if not isinstance(nb_refinements, int) or nb_refinements < 1:
            raise ValueError(f"Incorrect value for nb_refinements. ({nb_refinements})")

        if not isinstance(nb_polyorder, int) or nb_polyorder < 1:
            raise ValueError(f"Incorrect value for nb_polyorder. ({nb_polyorder})")

        if adaptivity not in {"disabled", "h", "p", "hp"}:
            raise ValueError(f"Incorrect value for adaptivity. ({adaptivity})")

        self.fields["electrostatic"]['analysis'] = analysis
        self.fields["electrostatic"]['solver'] = solver
        self.fields["electrostatic"]['matrix_solver'] = matrix_solver
        self.fields["electrostatic"]['nb_refinements'] = nb_refinements
        self.fields["electrostatic"]['nb_polyorder'] = nb_polyorder
        self.fields["electrostatic"]['adaptivity'] = adaptivity

    def add_field_magnetic(self, **kwargs):
        """
        :param analysis: 'steady', 'transient', 'harmonic'
        :param solver: 'linear', 'newton', 'picard'
        :param matrix_solver: 'umfpack', 'mumps', 'paraluation_it', 'paraluation_amg', 'mums_out'
        :param nb_refinements: Greater than 0 integer.
        :param nb_polyorder: Greater than 0 integer.
        :param adaptivity: 'disabled', 'h', 'p', 'hp'
        """

        analysis = kwargs.get("analysis", "steady")
        solver = kwargs.get("solver", "linear")
        matrix_solver = kwargs.get("matrix_solver", "umfpack")
        nb_refinements = kwargs.get("nb_refinements", 1)
        nb_polyorder = kwargs.get("nb_polyorder", 2)
        adaptivity = kwargs.get("adaptivity", "disabled")

        if analysis not in {"steady", "transient", "harmonic"}:
            raise ValueError(f"Magnetic field can can be: 'steady', 'transient', 'harmonic'. Got '{analysis}' instead.")

        if solver not in {"linear", "newton", "picard"}:
            raise ValueError(f"Solver type can be: 'linear', 'newton', 'picard'. Got {solver}")

        if matrix_solver not in {"umfpack", "mumps", "paraluation_it", "paraluation_amg", "mums_out"}:
            raise ValueError(f'Imporper matrix_solver type: "{matrix_solver}"')

        if not isinstance(nb_refinements, int) or nb_refinements < 1:
            raise ValueError(f"Incorrect value for nb_refinements. ({nb_refinements})")

        if not isinstance(nb_polyorder, int) or nb_polyorder < 1:
            raise ValueError(f"Incorrect value for nb_polyorder. ({nb_polyorder})")

        if adaptivity not in {"disabled", "h", "p", "hp"}:
            raise ValueError(f"Incorrect value for adaptivity. ({adaptivity})")

        self.fields["magnetic"]['analysis'] = analysis
        self.fields["magnetic"]['solver'] = solver
        self.fields["magnetic"]['matrix_solver'] = matrix_solver
        self.fields["magnetic"]['nb_refinements'] = nb_refinements
        self.fields["magnetic"]['nb_polyorder'] = nb_polyorder
        self.fields["magnetic"]['adaptivity'] = adaptivity

    def add_field_heat_transfer(self, **kwargs):
        """
        :param analysis: 'steady', 'transient'
        :param solver: 'linear', 'newton', 'picard'
        :param matrix_solver: 'umfpack', 'mumps', 'paraluation_it', 'paraluation_amg', 'mums_out'
        :param nb_refinements: Greater than 0 integer.
        :param nb_polyorder: Greater than 0 integer.
        :param adaptivity: 'disabled', 'h', 'p', 'hp'
        """

        analysis = kwargs.get("analysis", "steady")
        solver = kwargs.get("solver", "linear")
        matrix_solver = kwargs.get("matrix_solver", "umfpack")
        nb_refinements = kwargs.get("nb_refinements", 1)
        nb_polyorder = kwargs.get("nb_polyorder", 2)
        adaptivity = kwargs.get("adaptivity", "disabled")

        if analysis not in {"steady", "transient"}:
            raise ValueError(f"Heat transfer field can can be: 'steady', 'transient'. Got '{analysis}' instead.")

        if solver not in {"linear", "newton", "picard"}:
            raise ValueError(f"Solver type can be: 'linear', 'newton', 'picard'. Got {solver}")

        if matrix_solver not in {"umfpack", "mumps", "paraluation_it", "paraluation_amg", "mums_out"}:
            raise ValueError(f'Imporper matrix_solver type: "{matrix_solver}"')

        if not isinstance(nb_refinements, int) or nb_refinements < 1:
            raise ValueError(f"Incorrect value for nb_refinements. ({nb_refinements})")

        if not isinstance(nb_polyorder, int) or nb_polyorder < 1:
            raise ValueError(f"Incorrect value for nb_polyorder. ({nb_polyorder})")

        if adaptivity not in {"disabled", "h", "p", "hp"}:
            raise ValueError(f"Incorrect value for adaptivity. ({adaptivity})")

        self.fields["heat"]['analysis'] = analysis
        self.fields["heat"]['solver'] = solver
        self.fields["heat"]['matrix_solver'] = matrix_solver
        self.fields["heat"]['nb_refinements'] = nb_refinements
        self.fields["heat"]['nb_polyorder'] = nb_polyorder
        self.fields["heat"]['adaptivity'] = adaptivity

    def add_field_current_flow(self, **kwargs):
        """
        :param analysis: 'steady'
        :param solver: 'linear', 'newton'
        :param matrix_solver: 'umfpack', 'mumps', 'paraluation_it', 'paraluation_amg', 'mums_out'
        :param nb_refinements: Greater than 0 integer.
        :param nb_polyorder: Greater than 0 integer.
        :param adaptivity: 'disabled', 'h', 'p', 'hp'
        """

        analysis = kwargs.get("analysis", "steady")
        solver = kwargs.get("solver", "linear")
        matrix_solver = kwargs.get("matrix_solver", "umfpack")
        nb_refinements = kwargs.get("nb_refinements", 1)
        nb_polyorder = kwargs.get("nb_polyorder", 2)
        adaptivity = kwargs.get("adaptivity", "disabled")

        if analysis != "steady":
            raise ValueError(f"Current flow field can can be: 'steady'. Got '{analysis}' instead.")

        if solver not in {"linear", "newton"}:
            raise ValueError(f"Solver type can be: 'linear', 'newton'. Got {solver}")

        if matrix_solver not in {"umfpack", "mumps", "paraluation_it", "paraluation_amg", "mums_out"}:
            raise ValueError(f'Imporper matrix_solver type: "{matrix_solver}"')

        if not isinstance(nb_refinements, int) or nb_refinements < 1:
            raise ValueError(f"Incorrect value for nb_refinements. ({nb_refinements})")

        if not isinstance(nb_polyorder, int) or nb_polyorder < 1:
            raise ValueError(f"Incorrect value for nb_polyorder. ({nb_polyorder})")

        if adaptivity not in {"disabled", "h", "p", "hp"}:
            raise ValueError(f"Incorrect value for adaptivity. ({adaptivity})")

        self.fields["current"]['analysis'] = analysis
        self.fields["current"]['solver'] = solver
        self.fields["current"]['matrix_solver'] = matrix_solver
        self.fields["current"]['nb_refinements'] = nb_refinements
        self.fields["current"]['nb_polyorder'] = nb_polyorder
        self.fields["current"]['adaptivity'] = adaptivity

    def set_coordinate_type(self, coordinate_type):
        """
        This functions sets the type of the coordinate being used.

        :param coordinate_type: 'planar' or 'axisymmetric'
        """

        if coordinate_type.lower() in {"planar", "axisymmetric"}:
            self.coordinate_type = coordinate_type.lower()
        else:
            raise ValueError(f'There is no "{coordinate_type}" type of coordinate')

    def set_mesh_type(self, mesh_type):
        """
        This functions sets the type of the mesh being used.

        :param mesh_type: 'triangle', 'triangle_quad_fine_division', 'triangle_quad_rough_division',
        'triangle_quad_join', 'gmsh_triangle', 'gmsh_quad', 'gmsh_quad_delaunay'
        """

        if mesh_type.lower() in {
            "triangle",
            "triangle_quad_fine_division",
            "triangle_quad_rough_division",
            "triangle_quad_join",
            "gmsh_triangle",
            "gmsh_quad",
            "gmsh_quad_delaunay",
        }:
            self.mesh_type = mesh_type.lower()
        else:
            raise ValueError(f'There is no "{mesh_type}" type of mesh.')

    def add_boundary_condition(self, name, value, bctype, field):
        """
        Define a new boundary condition.

        :param name: name of the boundary condition
        :param value: The value of the boundary condition
        :param bctype: 'dirichlet', 'neumann'
        :param field: 'electrostatic', 'magnetic', 'heat', 'current'
        """
        if name not in self.boundary_conditions.keys():
            self.boundary_conditions[name] = (field, bctype, value)
        else:
            raise ValueError(f'There is a boundary condition called "{name}".')


if __name__ == "__main__":
    ag = Agros2DWrapper()
    ag.add_field_electrostatic()
    ag.add_boundary_condition("Source", 10, 'dirichlet', 'electrostatic')
    ag.add_boundary_condition("Ground", 0, 'dirichlet', 'electrostatic')
    ag.add_boundary_condition("Neumann", 0, 'neumann', 'electrostatic')
