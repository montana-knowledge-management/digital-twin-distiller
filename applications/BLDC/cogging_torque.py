from applications.BLDC.model import DIR_DATA
from model import *
import multiprocessing
from numpy import linspace
import matplotlib.pyplot as plt
from adze_modeler.utils import csv_write, csv_read


def calculate_cogging_torque():
    theta = linspace(0, 360 / 24, 501)
    models = [BLDCMotor(rotorangle=ti) for ti in theta]
    with multiprocessing.Pool(processes=4) as pool:
        T = pool.map(execute_model, models)
        csv_write(DIR_DATA / 'cogging_toruqe.csv', ['rotorangle', 'Torque'], theta, T)

def plot_cogging_torque():
    theta_sim, T_sim = csv_read(DIR_DATA / 'cogging_toruqe.csv')
    theta_ref, T_ref = csv_read(DIR_DATA / 'cogging_ref.csv')

    Tmax = max(T_sim)
    idx = T_sim.index(Tmax)
    print(f'Tmax: {Tmax:.3f} Nm @ {theta_sim[idx]} °')

    plt.figure()
    plt.plot(theta_sim, T_sim, 'b-', label='Simulation')
    plt.plot(theta_ref, T_ref, 'r-', label='Reference')
    
    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Rotor angle [°]")
    plt.ylabel("Coggin Torque [Nm]")
    plt.legend()
    plt.savefig(DIR_MEDIA / "cogging_toruqe.pdf", bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    from adze_modeler.utils import setup_matplotlib
    setup_matplotlib()

    # calculate_cogging_torque()
    plot_cogging_torque()
