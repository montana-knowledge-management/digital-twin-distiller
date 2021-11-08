from digital_twin_distiller.utils import getID


class BoundaryCondition:
    accepted_keys = {
        "electrostatic": [],
        "magnetic": [],
        "heat": [],
        "current": [],
    }

    def __init__(self, name, field_type):
        """
        :param name: name of the boundary condition
        :param field_type: 'electrostatic', 'magnetic', 'heat', 'current'
        """
        self.name = name
        self.field = field_type
        self.valuedict = {}
        self.assigned = set()
        self.type = None

        if self.field not in {"electrostatic", "magnetic", "heat", "current"}:
            raise ValueError(
                f'There is no "{field_type}" field type. Accepted values are '
                f'"electrostatic", "magnetic", "heat", "current".'
            )

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
        "heat": ["temperature"],
        "current": ["fixed_voltage"],
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
        "heat": [
            "heat_flux",
            "heat_transfer_coeff",
            "convection",
            "emissivity",
            "radiation",
        ],
        "current": ["current_density"],
    }

    def __init__(self, name, field_type, **kwargs):
        super().__init__(name, field_type)
        self.type = "neumann"

        for key, value in kwargs.items():
            self.set_value(key, value)


class PeriodicBoundaryCondition(BoundaryCondition):
    accepted_keys = {
        "electrostatic": [],
        "magnetic": [],
        "heat": [],
        "current": [],
    }

    def __init__(self, name, field_type, **kwargs):
        super().__init__(name, field_type)
        self.type = "neumann"


class AntiPeriodicBoundaryCondition(BoundaryCondition):
    accepted_keys = {
        "electrostatic": [],
        "magnetic": [],
        "heat": [],
        "current": [],
    }

    def __init__(self, name, field_type, **kwargs):
        super().__init__(name, field_type)
        self.type = "neumann"


class AntiPeriodicAirGap(BoundaryCondition):
    accepted_keys = {
        "electrostatic": [],
        "magnetic": [],
        "heat": [],
        "current": [],
    }

    def __init__(self, name, field_type, **kwargs):
        super().__init__(name, field_type)
        self.type = "antiperiodic"
        # backward compatibility
        if "inner_angle" in kwargs.keys():
            self.angle = kwargs["inner_angle"]
        elif "outer_angle" in kwargs.keys():
            self.angle = kwargs["outer_angle"]
        else:
            self.angle = kwargs.get("angle", 0)


class PeriodicAirGap(BoundaryCondition):
    accepted_keys = {
        "electrostatic": [],
        "magnetic": [],
        "heat": [],
        "current": [],
    }

    def __init__(self, name, field_type, **kwargs):
        super().__init__(name, field_type)
        self.type = "periodic"
        # backward compatibility
        if "inner_angle" in kwargs.keys():
            self.angle = kwargs["inner_angle"]
        elif "outer_angle" in kwargs.keys():
            self.angle = kwargs["outer_angle"]
        else:
            self.angle = kwargs.get("angle", 0)
