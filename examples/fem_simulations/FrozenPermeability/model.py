import operator as op
from copy import copy
from itertools import cycle
from math import acos, atan2, cos, degrees, fmod, hypot, pi, radians, sin

import numpy as np

from digital_twin_distiller import ModelPiece
from digital_twin_distiller.boundaries import DirichletBoundaryCondition
from digital_twin_distiller.material import Material
from digital_twin_distiller.metadata import FemmMetadata
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.objects import Node
from digital_twin_distiller.platforms.femm import Femm
from digital_twin_distiller.snapshot import Snapshot

ModelDir.set_base(__file__)


def cart2pol(x: float, y: float):
    r = hypot(x, y)
    phi = degrees(atan2(y, x))
    phi = fmod(phi + 1e-3, 360.0) - 1e-3
    return r, phi


def pol2cart(rho: float, phi: float):
    x = rho * cos(radians(phi))
    y = rho * sin(radians(phi))
    return x, y


def project(ax: float, ay: float, bx: float, by: float):
    """
    project vector 'a' onto vector 'b'.
    """
    c = ax * bx + ay * by / (bx**2 + by**2)
    return c * bx, c * by


class FrozenPermeability(BaseModel):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._init_directories()

        self.I0 = kwargs.get("I0", 0.0)
        self.alpha = kwargs.get("alpha", 0.0)
        self.rotorangle = kwargs.get("rotorangle", 0.0)

        self.Nturns = 152
        self.coil_area = 1
        self.J0 = self.I0 * self.Nturns / self.coil_area

        self.JU = self.J0 * cos(radians(self.alpha))
        self.JV = self.J0 * cos(radians(self.alpha + 120))
        self.JW = self.J0 * cos(radians(self.alpha + 240))

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False
        femm_metadata.depth = 1000  # TODO: Correct this value

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        # creating default materials
        airgap = Material("airgap")
        air = Material("air")
        steel_m19 = Material("M-19")
        steel_1018 = Material("1018")
        magnet = Material("Magnet")
        coil = Material("coil")

        # default material configuration
        airgap.meshsize = 0.1

        steel_m19.thickness = 0.635
        steel_m19.lamination_type = "inplane"
        steel_m19.fill_factor = 0.98
        steel_m19.conductivity = 1.9e6
        steel_m19.b = [
            0.0,
            0.05,
            0.10,
            0.15,
            0.20,
            0.25,
            0.30,
            0.35,
            0.40,
            0.45,
            0.50,
            0.55,
            0.60,
            0.65,
            0.70,
            0.75,
            0.80,
            0.85,
            0.90,
            0.95,
            1.00,
            1.05,
            1.10,
            1.15,
            1.20,
            1.25,
            1.30,
            1.35,
            1.40,
            1.45,
            1.50,
            1.55,
            1.60,
            1.65,
            1.70,
            1.75,
            1.80,
            1.85,
            1.90,
            1.95,
            2.00,
            2.05,
            2.10,
            2.15,
            2.20,
            2.25,
            2.30,
        ]
        steel_m19.h = [
            0.0,
            15.120714,
            22.718292,
            27.842733,
            31.871434,
            35.365044,
            38.600588,
            41.736202,
            44.873979,
            48.087807,
            51.437236,
            54.975221,
            58.752993,
            62.823644,
            67.245285,
            72.084406,
            77.420100,
            83.350021,
            89.999612,
            97.537353,
            106.201406,
            116.348464,
            128.547329,
            143.765431,
            163.754169,
            191.868158,
            234.833507,
            306.509769,
            435.255202,
            674.911968,
            1108.325569,
            1813.085468,
            2801.217421,
            4053.653117,
            5591.106890,
            7448.318413,
            9708.815670,
            12486.931615,
            16041.483644,
            21249.420624,
            31313.495878,
            53589.446877,
            88477.484601,
            124329.410540,
            159968.569300,
            197751.604272,
            234024.751347,
        ]

        steel_1018.conductivity = 5.8e6
        steel_1018.b = [0.0, 0.2503, 0.925, 1.250, 1.390, 1.525, 1.71, 1.87, 1.955, 2.02, 2.11, 2.225, 2.43]
        steel_1018.h = [
            0.0,
            238.7325,
            795.775,
            1591.55,
            2387.325,
            3978.875,
            7957.75,
            15915.5,
            23873.25,
            39788.75,
            79577.5,
            159155,
            318310,
        ]
        steel_1018.phi_hmax = 20

        magnet.mu_r = 1.04496
        magnet.coercivity = 891000

        coil.diameter = 0.643808017739015
        coil.lamination_type = "magnetwire"
        coil.conductivity = 58e6

        # Generating additional materials
        # with the above base materials

        # Coils
        # PHASE U
        phase_U_positive = copy(coil)
        phase_U_positive.name = "U+"
        phase_U_positive.Je = self.JU

        phase_U_negative = copy(coil)
        phase_U_negative.name = "U-"
        phase_U_negative.Je = -self.JU

        # PHASE V
        phase_V_positive = copy(coil)
        phase_V_positive.name = "V+"
        phase_V_positive.Je = self.JV

        phase_V_negative = copy(coil)
        phase_V_negative.name = "V-"
        phase_V_negative.Je = -self.JV

        # PHASE W
        phase_W_positive = copy(coil)
        phase_W_positive.name = "W+"
        phase_W_positive.Je = self.JW

        phase_W_negative = copy(coil)
        phase_W_negative.name = "W-"
        phase_W_negative.Je = -self.JW

        # magnets
        flipper = cycle((0, 180))
        for i in range(6):
            mi = copy(magnet)
            mi.name = f"{mi.name}_{i+1}"
            mi.remanence_angle = fmod(180 + i * 360 / 6 + next(flipper) + self.rotorangle, 360.0)
            self.snapshot.add_material(mi)
            label_x, label_y = pol2cart(28.5, i * 360 / 6 + self.rotorangle)
            self.snapshot.assign_material(label_x, label_y, mi.name)

        # coil labels
        winding_labels = cycle(("U+", "W-", "V+", "U-", "W+", "V-"))
        for i in range(18):
            label_x, label_y = pol2cart(40, i * 360 / 18)
            self.assign_material(label_x, label_y, next(winding_labels))

        self.assign_material(0, 30.25 + 0.75 / 2, "airgap")

        self.snapshot.add_material(air)
        self.snapshot.add_material(airgap)
        self.snapshot.add_material(steel_m19)
        self.snapshot.add_material(steel_1018)
        self.snapshot.add_material(phase_U_positive)
        self.snapshot.add_material(phase_U_negative)
        self.snapshot.add_material(phase_V_positive)
        self.snapshot.add_material(phase_V_negative)
        self.snapshot.add_material(phase_W_positive)
        self.snapshot.add_material(phase_W_negative)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

        self.assign_boundary_arc(0, 53, "a0")
        self.assign_boundary_arc(0, -53, "a0")

    def add_postprocessing(self):
        r = 30.25 + 0.75 / 2
        phi = np.linspace(0, 2 * np.pi, 1001)
        for phi_i in phi:
            self.snapshot.add_postprocessing("point_value", (r * np.cos(phi_i), r * np.sin(phi_i)), "Bx")
            self.snapshot.add_postprocessing("point_value", (r * np.cos(phi_i), r * np.sin(phi_i)), "By")

    def build_geometry(self):
        stator = ModelPiece("stator")
        stator.load_piece_from_dxf(ModelDir.RESOURCES / "stator.dxf")

        rotor = ModelPiece("rotor")
        rotor.load_piece_from_dxf(ModelDir.RESOURCES / "rotor.dxf")
        rotor.rotate(alpha=self.rotorangle)

        self.geom.merge_geometry(stator.geom)
        self.geom.merge_geometry(rotor.geom)

        for i in range(len(self.geom.circle_arcs)):
            self.geom.circle_arcs[i].max_seg_deg = 1

        self.snapshot.add_geometry(self.geom)

        self.assign_material(0, 0, "air")
        self.assign_material(0, 50, "M-19")
        self.assign_material(0, 20, "1018")

    def __call__(self, cleanup=True, devmode=False, timeout=10000):
        r_ = super().__call__(cleanup=cleanup, devmode=devmode, timeout=timeout)

        Bx = list(map(op.itemgetter(2), r_["Bx"]))
        By = list(map(op.itemgetter(2), r_["By"]))
        x = list(map(op.itemgetter(0), r_["Bx"]))
        y = list(map(op.itemgetter(1), r_["By"]))

        res = {"Br": [], "Bphi": [], "r": [], "phi": []}
        for Bxi, Byi, xi, yi in zip(Bx, By, x, y):
            r_i, phi_i = cart2pol(xi, yi)

            et = Node(xi, yi).rotate(radians(90))
            et = et / abs(et)
            er = Node.from_polar(1, phi_i)

            # Itt kell kiszámítani a B(Bx, By) vektor vetületét a tangenciális (et) meg a radiális (er) irányú
            # egységvektorokra
            B_ri = project(Bxi, Byi, er.x, er.y)
            B_ri = hypot(*B_ri)

            B_phi = project(Bxi, Byi, et.x, et.y)
            B_phi = hypot(*B_phi)

            res["r"].append(r_i)
            res["phi"].append(phi_i)
            res["Br"].append(B_ri)
            res["Bphi"].append(B_phi)

        res["x"] = x
        res["y"] = y
        res["Bx"] = Bx
        res["By"] = By

        return res


if __name__ == "__main__":
    import json

    import matplotlib.pyplot as plt

    m = FrozenPermeability(I0=0.0, alpha=0.0, rotorangle=0.0, exportname="dev")
    res_ = m(cleanup=False, devmode=False)

    with open(ModelDir.DATA / "example.json", "w", encoding="utf-8") as f:
        json.dump(res_, f, ensure_ascii=True, indent=2)

    with open(ModelDir.DATA / "example.json", encoding="utf-8") as f:
        res_ = json.load(f)

    plt.figure(figsize=(10, 10))
    plt.plot(res_["phi"][1:], res_["Bx"][1:])
    plt.plot(res_["phi"][1:], res_["By"][1:])

    plt.savefig(ModelDir.MEDIA / "example.png", bbox_inches="tight", dpi=450)
    plt.show()
