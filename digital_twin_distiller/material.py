from copy import deepcopy


class Material:
    def __init__(self, name):
        self.name = name
        self.mu_r = 1.0
        self.epsioln_r = 1.0
        self.conductivity = 0.0
        self.b = []
        self.h = []
        self.Je = 0.0  # External current density, can be complex
        self.Rho = 0.0  # Volume charge density
        self.remanence_angle = 0.0
        self.remanence = 0.0
        self.coercivity = 0.0
        self.angluar_velocity = 0.0
        self.vx = 0.0
        self.vy = 0.0

        # Femm realted
        self.thickness = 0
        self.lamination_type = 0
        self.fill_factor = 0
        self.diameter = 1.0
        self.phi_hmax = 0.0

        self.assigned = []  # a set of (x, y) tuples
        self.meshsize = 0

    def __copy__(self):
        # newmaterial = Material(self.name)
        # newmaterial.name = self.name
        # newmaterial.mu_r = self.mu_r
        # newmaterial.epsioln_r = self.epsioln_r
        # newmaterial.conductivity = self.conductivity
        # newmaterial.b = self.b.copy() or []
        # newmaterial.h = self.h.copy() or []
        # newmaterial.Je = self.Je
        # newmaterial.Rho = self.Rho
        # return newmaterial

        return deepcopy(self)
