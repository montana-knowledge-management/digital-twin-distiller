import csv
import math
import multiprocessing
import random
import uuid
from pathlib import Path
from statistics import fmean

import matplotlib.pyplot as plt
from model import PriusMotor
from numpy.core.function_base import linspace


class ParametricPriusMotor(PriusMotor):
    def __init__(self, X: list, rotorangle: float = 0, exportname=None):
        super().__init__(rotorangle=rotorangle, exportname=exportname)
        assert len(X) == 5

        self.size_airgap = X[0]
        self.size_steel = X[1]
        self.size_coil = X[2]
        self.size_magnet = X[3]
        self.size_air = X[4]

    def define_materials(self):
        super().define_materials()
        self.snapshot.materials["air"].meshsize = self.size_air
        self.snapshot.materials["airgap"].meshsize = self.size_airgap
        self.snapshot.materials["U+"].meshsize = self.size_coil
        self.snapshot.materials["U-"].meshsize = self.size_coil
        self.snapshot.materials["V+"].meshsize = self.size_coil
        self.snapshot.materials["V-"].meshsize = self.size_coil
        self.snapshot.materials["W+"].meshsize = self.size_coil
        self.snapshot.materials["W-"].meshsize = self.size_coil
        self.snapshot.materials["M19_29G"].meshsize = self.size_steel
        self.snapshot.materials["M19_29GSF094"].meshsize = self.size_steel
        self.snapshot.materials["N36Z_50_r"].meshsize = self.size_magnet
        self.snapshot.materials["N36Z_50_l"].meshsize = self.size_magnet

    def add_postprocessing(self):
        super().add_postprocessing()
        self.snapshot.add_postprocessing("mesh_info", None, None)

    def __repr__(self):
        return f"{self.rotorangle:.3f} {self.size_airgap} {self.size_steel} {self.size_coil} {self.size_magnet} {self.size_air}"


def execute_model(model: ParametricPriusMotor):
    res = model(timeout=2000)
    torque = res["Torque"]
    nb_elements = int(res["elements"])
    print(f"{model} {torque:.3e}: - {nb_elements}")
    return model.rotorangle, torque, nb_elements


def analyze_cogging(X):
    dir_data = Path(__file__).parent / "data"
    rotorangles = linspace(0, 1, 51) * 45 / 8 * 2
    models = [ParametricPriusMotor(X, rotorangle=ri) for ri in rotorangles]

    res = []
    with multiprocessing.Pool(processes=4) as pool:
        res = pool.map(execute_model, models)

    nb_elements = int(fmean({ri[2] for ri in res}))
    max_T = max(res, key=lambda ri: ri[1])[1]

    filename = dir_data / "mesh_sensitivity" / f"{uuid.uuid4()}.csv"

    with open(filename, "w") as f:
        w = csv.writer(f)
        w.writerows(res)

    newcfg = {
        "size_airgap": X[0],
        "size_steel": X[1],
        "size_coil": X[2],
        "size_magnet": X[3],
        "size_air": X[4],
        "nb_elements": nb_elements,
        "peak_cogging": max_T,
        "filename": str(filename.relative_to(Path(__file__).parent)),
    }

    with open(dir_data / "mesh_sensitivity" / "config.csv", mode="a+") as f:
        writer = csv.DictWriter(f, fieldnames=newcfg.keys(), quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(newcfg)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    dir_data = Path(__file__).parent

    # get_x = lambda: [random.uniform(0.2, 1.5) for _ in range(5)]
    # analyze_cogging(get_x())

    cfglist = []
    with open(dir_data / "data" / "mesh_sensitivity" / "config.csv") as f:
        cfglist = list(csv.DictReader(f, quoting=csv.QUOTE_NONNUMERIC))

    cfglist = sorted(cfglist, key=lambda cfg_i: cfg_i["nb_elements"])
    NUM_COLORS = len(cfglist)

    nb_elements = []
    peak_cogging = []

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cm = plt.get_cmap("brg")
    ax.set_prop_cycle(color=[cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])

    for cfg_i in cfglist:
        theta = []
        T = []
        nb_elements.append(cfg_i["nb_elements"])
        peak_cogging.append(cfg_i["peak_cogging"] * 8)

        with open(dir_data / cfg_i["filename"]) as f:
            for theta_i, Ti, nb_elements_i in csv.reader(f, quoting=csv.QUOTE_NONNUMERIC):
                theta.append(theta_i)
                T.append(Ti * 8)

        plt.plot(theta, T, "-o", label=int(cfg_i["nb_elements"]))

    plt.grid()
    plt.xlabel("Rotor Angle [°]")
    plt.ylabel("Cogging Torque [Nm]")
    plt.legend()
    plt.show()

    plt.figure()
    plt.scatter(nb_elements, peak_cogging, c="b", alpha=0.7)
    plt.grid()
    plt.xlabel("Number of elements")
    plt.ylabel("Peak torque [Nm]")
    plt.show()
