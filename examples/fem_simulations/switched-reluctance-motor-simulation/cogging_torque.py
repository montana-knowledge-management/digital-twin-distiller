import multiprocessing

import matplotlib.pyplot as plt
from model import SRM, execute_model
from numpy import linspace

from digital_twin_distiller.modelpaths import *
from digital_twin_distiller.utils import csv_read, csv_write

ModelDir.set_base(__file__)


def calculate_cogging_torque():
    theta = linspace(0, 360 / 3, 80)
    models = [SRM(rotorangle=ti) for ti in theta]
    with multiprocessing.Pool(processes=4) as pool:
        T = pool.map(execute_model, models)
        csv_write(ModelDir.DATA / "cogging_torque.csv", ["rotorangle", "Torque"], theta, T)


def plot_cogging_torque():
    theta_sim, T_sim = csv_read(ModelDir.DATA / "cogging_torque.csv")

    Tmax = max(T_sim)
    idx = T_sim.index(Tmax)
    print(f"Tmax: {Tmax:.3f} Nm @ {theta_sim[idx]} °")

    plt.figure()
    plt.plot(theta_sim, T_sim, "b-", label="Simulation")

    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Rotor angle [°]")
    plt.ylabel("Cogging Torque [Nm]")
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "ct_validation.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":

    calculate_cogging_torque()
    plot_cogging_torque()
