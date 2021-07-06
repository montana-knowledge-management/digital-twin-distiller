from abc import ABCMeta
from abc import abstractmethod
from collections.abc import Iterable
from types import SimpleNamespace


def newline(n=1):
    return print("\n" * (n - 1))


class Field(metaclass=ABCMeta):
    def __init__(self):
        self.name = None
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

    def export(self):
        print(f"# {self.name}")
        print(f'{self.name} = a2d.field("{self.name}")')
        print(f'{self.name}.analysis_type = "{self.analysis}"')
        print(f"{self.name}.number_of_refinements = {self.nb_refinements}")
        print(f"{self.name}.polynomial_order = {self.polyorder}")
        print(f'{self.name}.adaptivity_type = "{self.adaptivity}"')
        print(f'{self.name}.solver = "{self.solver}"')

        newline(2)
        print("# boundaries")
        self._export_boundary_conditions()

        newline(2)
        print("# materials")
        self._export_materials()

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
    def _export_materials(self):
        ...

    @abstractmethod
    def _export_boundary_conditions(self):
        ...


class ElectrostaticField(Field):
    def __init__(self):
        super().__init__()
        self.name = "electrostatic"

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

    def _export_materials(self):
        for name, mi in self.materials.items():
            print(
                f'{self.name}.add_material("{name}", '
                f'{{"electrostatic_permittivity" : {mi[0]}, '
                f'"electrostatic_charge_density" : {mi[1]}}})'
            )

    def _export_boundary_conditions(self):
        for key, bi in self.bc_dirichlet.items():
            print(
                f'{self.name}.add_boundary("{key}", '
                f'"electrostatic_potential", '
                f'{{"electrostatic_potential" : {bi}}})'
            )

        for key, bi in self.bc_neumann.items():
            print(
                f'{self.name}.add_boundary("{key}", '
                f'"electrostatic_surface_charge_density", '
                f'{{"electrostatic_surface_charge_density" : {bi}}})'
            )


class MagneticField(Field):
    def __init__(self):
        super().__init__()
        self.name = "magnetic"

    def set_analysis(self, analysis_type):
        if analysis_type == "steady":
            self.analysis = "steadystate"

        elif analysis_type == "transient":
            self.analysis = "transient"

        elif analysis_type == "harmonic":
            self.analysis = "harmonic"

        else:
            raise ValueError(f'Accepted values: "steady", "transient", "harmonic". Got "{analysis_type}"')

    def set_solver(self, solver_type):

        if solver_type == "linear":
            self.solver = "linear"
        elif solver_type == "picard":
            self.solver = "picard"
        elif solver_type == "newton":
            self.solver = "newton"
        else:
            raise ValueError(f'Accepted values: "linear", "picard", "newton". Got "{solver_type}"')

    def add_boundary_condition(self, type_, name, value=0):
        if type_ not in {"d", "n"}:
            raise ValueError(f"Expected values for type_ are 'd' (dirichlet) or 'n' (neumann). Got '{type_}'.")

        if type_ == "d":
            self.bc_dirichlet[name] = value
        else:
            self.bc_neumann[name] = value

    def add_material(
        self,
        name,
        current_density_external=0,
        total_current=0,
        conductivity=0,
        velocity_x=0,
        velocity_y=0,
        velocity_angular=0,
        mur=1,
        remanence_angle=0,
        remanence=0,
        H=None,
        B=None,
        **kwargs,
    ):
        mu0 = 4 * (1146408 / 364913) * 10 ** -7
        newmaterial = SimpleNamespace(
            Ipresc=False,
            Je=current_density_external,
            Itot=total_current,
            Sigma=conductivity,
            vx=velocity_x,
            vy=velocity_y,
            omega=velocity_angular,
            mur=mur,
            B=None,
            remanence=remanence,
            remanence_angle=remanence_angle,
        )

        if total_current:
            newmaterial.total_current_prescribed = True

        if B and H and isinstance(B, Iterable) and isinstance(H, Iterable) and len(B) == len(H):
            self.set_solver("picard")
            mu_r = list()
            for Bi, Hi in zip(B, H):
                mu_r.append(Bi / (Hi * mu0))

            interpolation_method = kwargs.get("interpolation", "piecewise")
            newmaterial.interpolation = interpolation_method
            newmaterial.mur = mu_r.copy()
            newmaterial.B = B.copy()

        self.materials[name] = newmaterial

    def _export_materials(self):
        # TODO: check proper syntax for nonlinear materials
        for name, mi in self.materials.items():
            mdict = {
                "magnetic_remanence_angle": mi.remanence_angle,
                "magnetic_velocity_y": mi.vy,
                "magnetic_current_density_external_real": mi.Je.real,
                "magnetic_current_density_external_imag": mi.Je.imag,
                "magnetic_permeability": mi.mur,
                "magnetic_conductivity": mi.Sigma,
                "magnetic_total_current_prescribed": mi.Ipresc,
                "magnetic_remanence": mi.remanence,
                "magnetic_velocity_angular": mi.omega,
                "magnetic_total_current_real": mi.Itot.real,
                "magnetic_velocity_x": mi.vx,
                "magnetic_total_current_imag": mi.Itot.imag,
            }
            if self.analysis != "harmonic":
                mdict.pop("magnetic_total_current_imag")
                mdict.pop("magnetic_current_density_external_imag")

            print(f'magnetic.add_material("{name}", {str(mdict)})')

    def _export_boundary_conditions(self):
        for name, bi in self.bc_dirichlet.items():
            bcdict = {
                "magnetic_potential_real": bi.real,
                # "magnetic_surface_current_real": 0,
            }
            if self.analysis == "harmonic":
                bcdict["magnetic_potential_imag"] = bi.imag

            print(f'{self.name}.add_boundary("{name}", "magnetic_potential", {str(bcdict)})')

            for name, bi in self.bc_neumann.items():
                bcdict = {"magnetic_surface_current_real": bi.real}
                if self.analysis == "harmonic":
                    bcdict["magnetic_surface_current_imag"] = bi.imag

                print(f'{self.name}.add_boundary("{name}", "magnetic_surface_current", {str(bcdict)})')


class HeatFlowField(Field):
    def __init__(self):
        super().__init__()
        self.name = "heat"
