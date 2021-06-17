from abc import ABCMeta, abstractmethod


def newline(n):
    return print("\n" * (n - 1))


class Field(metaclass=ABCMeta):
    def __init__(self):
        self.analysis = "steadystate"
        self.solver = "linear"
        self.matrix_solver = "umfpack"
        self.nb_refinements = 1
        self.polyorder = 2
        self.adaptivity = "disabled"
        self.bc_dirichlet = {}
        self.bc_neumann = {}
        self.materials = {}

    def set_refinements(self, n):
        if isinstance(n, int) and n > 0:
            self.nb_refinements = n
        else:
            raise ValueError(f"Number of refinement is a positive integer. Got {n}.")

    def set_polynomial_order(self, n):
        if isinstance(n, int) and n > 0:
            self.polyorder = n
        else:
            raise ValueError(f"Polynomial order is a positive integer. Got {n}.")

    def set_adaptivity(self, adaptivity_type="disabled"):
        if adaptivity_type == "disabled":
            self.adaptivity = "disabled"
        elif adaptivity_type == "h":
            self.adaptivity = "h"
        elif adaptivity_type == "p":
            self.adaptivity = "p"
        elif adaptivity_type == "hp":
            self.adaptivity = "hp"
        else:
            raise ValueError(
                f'Accepted values for adaptivity_type: "disabled", "h", "p", "hp". Got "{adaptivity_type}"'
            )

    @abstractmethod
    def set_analysis(self, analysis_type):
        ...

    @abstractmethod
    def set_solver(self, solver_type):
        ...

    @abstractmethod
    def add_boundary_condition(self, name):
        ...

    @abstractmethod
    def add_material(self, name):
        ...

    @abstractmethod
    def export(self):
        ...


class ElectrostaticField(Field):
    def __init__(self):
        super().__init__()

    def set_analysis(self, analysis_type):
        if analysis_type in {"steady", "steadystate", "stationary", "static"}:
            self.analysis = "steadystate"
        else:
            raise ValueError(f'Electrostatic problems can be "steady" only. Got "{analysis_type}".')

    def set_solver(self, solver_type):
        if solver_type == "linear":
            self.solver = "linear"
        else:
            raise ValueError(f'Electrostatic solver is "linear" only. Got "{solver_type}".')

    def add_boundary_condition(self, type_, name, value=0):
        """
        :param type_: 'd' - dirichlet, 'n' - neumann
        :param name: name of the boundary condition
        :param value: the value of the boundary condition
        """
        if type_ not in {"d", "n"}:
            raise ValueError(f"Expected values for type_ are 'd' (dirichlet) or 'n' (neumann). Got '{type_}'.")

        if type_ == "d":
            self.bc_dirichlet[name] = value
        else:
            self.bc_neumann[name] = value

    def add_material(self, name, epsilon_r=1, sigma=0):
        self.materials[name] = (epsilon_r, sigma)

    def export(self):
        print("# Electrostatic")
        print(f'electrostatic = a2d.field("electrostatic")')
        print(f'electrostatic.analysis_type = "{self.analysis}"')
        print(f"electrostatic.number_of_refinements = {self.nb_refinements}")
        print(f"electrostatic.polynomial_order = {self.polyorder}")
        print(f'electrostatic.adaptivity_type = "{self.adaptivity}"')
        print(f'electrostatic.solver = "{self.solver}"')
        newline(2)
        print("# boundaries")
        for key, bi in self.bc_dirichlet.items():
            print(
                f'electrostatic.add_boundary("{key}", '
                f'"electrostatic_potential", '
                f'{{"electrostatic_potential" : {bi}}})'
            )

        for key, bi in self.bc_neumann.items():
            print(
                f'electrostatic.add_boundary("{key}", '
                f'"electrostatic_surface_charge_density", '
                f'{{"electrostatic_surface_charge_density" : {bi}}})'
            )

        newline(2)
        print("# materials")
        for name, mi in self.materials.items():
            print(f'electrostatic.add_material("{name}", '
                  f'{{"electrostatic_permittivity" : {mi[0]}, '
                  f'"electrostatic_charge_density" : {mi[1]}}})')