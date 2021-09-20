import csv
import math
import multiprocessing
import operator
import random
import uuid
from pathlib import Path
from statistics import fmean

import matplotlib.pyplot as plt
from model import PriusMotor
from numpy import poly1d
from numpy import polyfit
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
        self.snapshot.materials["M19_29GSF094"].meshsize = self.size_steel
        self.snapshot.materials["N36Z_50_r"].meshsize = self.size_magnet
        self.snapshot.materials["N36Z_50_l"].meshsize = self.size_magnet

    def add_postprocessing(self):
        super().add_postprocessing()
        self.snapshot.add_postprocessing("mesh_info", None, None)

    def __repr__(self):
        return f"{self.rotorangle:.3f}° {self.size_airgap} {self.size_steel} {self.size_coil} {self.size_magnet} {self.size_air}"


def execute_model(model: ParametricPriusMotor):
    res = model(timeout=2000)
    torque = res["Torque"]
    nb_elements = int(res["elements"])
    print(f"{model}: {torque*8:.3f} Nm - {nb_elements}")
    return model.rotorangle, torque, nb_elements


def analyze_cogging(X):
    X = [round(float(xi), 3) for xi in X]
    print(X)
    dir_data = Path(__file__).parent / "data"
    rotorangles = linspace(0, 360 / 48, 31)
    models = [ParametricPriusMotor(X, rotorangle=ri) for ri in rotorangles]

    res = []
    with multiprocessing.Pool(processes=4) as pool:
        res = pool.map(execute_model, models)

    nb_elements = int(fmean({ri[2] for ri in res}))
    max_T = max(res, key=lambda ri: ri[1])[1]
    min_T = min(res, key=lambda ri: ri[1])[1]
    T_pp = max_T * 8 - min_T * 8

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
        "peak_cogging": T_pp,
        "filename": str(filename.relative_to(Path(__file__).parent)),
    }

    with open(dir_data / "mesh_sensitivity" / "config.csv", mode="a+") as f:
        writer = csv.DictWriter(f, fieldnames=newcfg.keys(), quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(newcfg)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    dir_data = Path(__file__).parent

    # get_x = lambda: [random.uniform(0.1, 0.5) for _ in range(5)]
    # analyze_cogging(get_x())
    # analyze_cogging([1, 0.1, 1, 1, 1])

    #
    # analyze_cogging([1, 1, 1, 1, 1])
    # analyze_cogging([0.6, 1, 1, 1, 1])
    # analyze_cogging([0.5, 1, 1, 1, 1])
    # analyze_cogging([0.4, 1, 1, 1, 1])
    # analyze_cogging([0.3, 1, 1, 1, 1])
    # analyze_cogging([0.2, 1, 1, 1, 1])
    # analyze_cogging([0.1, 1, 1, 1, 1])

    cfglist = []
    with open(dir_data / "data" / "mesh_sensitivity" / "config.csv") as f:
        cfglist = list(csv.DictReader(f, quoting=csv.QUOTE_NONNUMERIC))

    cfglist = sorted(cfglist, key=lambda cfg_i: cfg_i["nb_elements"])
    NUM_COLORS = len(cfglist)

    nb_elements = []
    rms_torque = []

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cm = plt.get_cmap("binary")
    ax.set_prop_cycle(color=[cm(1.0 * i / NUM_COLORS) for i in range(2, NUM_COLORS + 2)])

    N = len(cfglist)
    for i, cfg_i in enumerate(cfglist):
        theta = []
        T = []
        nb_elements.append(cfg_i["nb_elements"])

        with open(dir_data / cfg_i["filename"]) as f:
            for theta_i, Ti, nb_elements_i in csv.reader(f, quoting=csv.QUOTE_NONNUMERIC):
                theta.append(theta_i)
                T.append(Ti * 8)

        T_rms = math.sqrt(sum(map(lambda ti: ti ** 2, T)) / len(T))
        rms_torque.append(T_rms)

        p = poly1d(polyfit(theta, T, 11))
        theta_fine = linspace(min(theta), max(theta), 501)

        if i == 0:
            plt.plot(theta_fine, p(theta_fine), "b-", label=int(cfg_i["nb_elements"]), linewidth=2)
            # plt.plot(theta, T, "r-", label=int(cfg_i["nb_elements"]))
        elif i == N - 1:
            plt.plot(theta_fine, p(theta_fine), "r-", label=int(cfg_i["nb_elements"]), linewidth=2)
        else:
            plt.plot(theta_fine, p(theta_fine), linestyle="dashdot", label=int(cfg_i["nb_elements"]))

    plt.grid()
    plt.xlabel("Rotor Angle [°]")
    plt.ylabel("Cogging Torque [Nm]")
    plt.legend()
    plt.show()

    plt.figure()
    plt.scatter(nb_elements, rms_torque, c="b", alpha=0.7)
    plt.grid()
    plt.xlabel("Number of elements")
    plt.ylabel("RMS torque [Nm]")
    plt.show()
