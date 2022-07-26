from copy import deepcopy


class Material:
    def __init__(self, name, **kwargs):
        self.name = name
        self.mu_r = kwargs.get("mu_r", 1.0)
        self.epsioln_r = kwargs.get("epsioln_r", 1.0)
        self.conductivity = kwargs.get("conductivity", 0.0)
        self.b = kwargs.get("b", []).copy()
        self.h = kwargs.get("h", []).copy()
        self.Je = kwargs.get("Je", 0.0)  # External current density, can be complex
        self.Rho = kwargs.get("Rho", 0.0)  # Volume charge density
        self.remanence_angle = kwargs.get("remanence_angle", 0.0)
        self.remanence = kwargs.get("remanence", 0.0)
        self.coercivity = kwargs.get("coercivity", 0.0)
        self.angluar_velocity = kwargs.get("angluar_velocity", 0.0)
        self.vx = kwargs.get("vx", 0.0)
        self.vy = kwargs.get("vy", 0.0)

        # Femm related
        self.thickness = kwargs.get("thickness", 0)
        self.lamination_type = kwargs.get("lamination_type", 0)
        self.fill_factor = kwargs.get("fill_factor", 0)
        self.diameter = kwargs.get("diameter", 1.0)
        self.phi_hmax = kwargs.get("phi_hmax", 0.0)

        # FEMM HEAT
        self.kx = kwargs.get("kx", 1.0)
        self.ky = kwargs.get("ky", 1.0)
        self.qv = kwargs.get("qv", 0.0)
        self.kt = kwargs.get("kt", 0.0)

        # AGROS2D related
        # HEAT
        self.material_density = kwargs.get("material_density", 0.0)
        self.heat_conductivity = kwargs.get("heat_conductivity", 385.0)
        self.volume_heat = kwargs.get("volume_heat", 0.0)
        self.specific_heat = kwargs.get("specific_heat", 0.0)

        self.assigned = []  # a set of (x, y) tuples
        self.meshsize = kwargs.get("meshsize", 0)

    def __copy__(self):
        return deepcopy(self)
