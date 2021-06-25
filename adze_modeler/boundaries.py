from adze_modeler.utils import getID


class BoundaryCondition:
    accepted_keys = {
        "electrostatic": [],
        "magnetic": [],
        "heat": [],
        "current": [],
    }

    def __init__(self, name, field_type, **kwargs):
        """
        :param name: name of the boundary condition
        :param field_type: 'electrostatic', 'magnetic', 'heat', 'current'
        """
        self.name = name
        self.field = field_type
        self.valuedict = {}
        self.assigned = set()
        self.type = None

        # setting initial values to valuedict
        # for key in self.accepted_keys[self.field]:
        #     self.valuedict[key] = 0.0

    def set_value(self, key, value):
        if key in self.accepted_keys[self.field]:
            self.valuedict[key] = value
        else:
            raise ValueError(f'There is no "{key}" in {self.field} dirichlet boundary condition.')

    def __str__(self):
        st = f"name: {self.name}, type: {self.field}-{self.type}, value(s): "
        for key in self.valuedict.keys():
            st += f"{key}: {self.valuedict[key]}"
        return st


class DirichletBoundaryCondition(BoundaryCondition):
    accepted_keys = {
        "electrostatic": ["fixed_voltage"],
        "magnetic": ["magnetic_potential"],
        "heat": [],
        "current": [],
    }

    def __init__(self, name, field_type, **kwargs):
        super().__init__(name, field_type)
        self.type = "dirichlet"

        for key, value in kwargs.items():
            self.set_value(key, value)


class NeumannBoundaryCondition(BoundaryCondition):
    accepted_keys = {
        "electrostatic": ["surface_charge_density"],
        "magnetic": ["surface_current"],
        "heat": [],
        "current": [],
    }

    def __init__(self, name, field_type, **kwargs):
        super().__init__(name, field_type)
        self.type = "neumann"
