import multiprocessing
from pathlib import Path

import matplotlib.pyplot as plt
from model import PriusMotor
from numpy import linspace
from cycler import cycler

"""
http://phdengineeringem.blogspot.com/2018/06/toyota-prius-femm-cogging-torque.html

https://magneticsmag.com/simulating-the-toyota-prius-electric-motor/

https://au-193-34.rev.ensam.fr/bitstream/handle/10985/6752/JEPE11090903-ES-final.pdf?sequence=1&isAllowed=y

"""


def execute_model(model_i: PriusMotor):
    res = model_i(cleanup=False)
    T = res["Torque"]
    print(f"{abs(model_i.rotorangle):.2f} ° - {T*8:.3f}")
    return T


def calculate_cogging_torque():
    theta = linspace(0, 360 / 48, 51)
    models = [PriusMotor(rotorangle=-theta_i) for theta_i in theta]
    # models = [PriusMotor(alpha=-theta_i, I0=I0) for theta_i in theta]

    with multiprocessing.Pool(processes=4) as pool:
        res = pool.map(execute_model, models)

    with open(models[0].dir_data / f"cogging_torque.csv", "w") as f:
        for theta_i, T_i in zip(theta, res):
            f.write(f"{theta_i}, {T_i}\n")


def calculate_torque(I0=0):
    theta = linspace(0, 180, 31)
    # models = [PriusMotor(rotorangle=theta_i, I0=250) for theta_i in theta]
    models = [PriusMotor(alpha=-theta_i, I0=I0) for theta_i in theta]

    with multiprocessing.Pool(processes=4) as pool:
        res = pool.map(execute_model, models)

    with open(models[0].dir_data / "torque" / f"torque_{int(I0)}.csv", "w") as f:
        for theta_i, T_i in zip(theta, res):
            f.write(f"{theta_i}, {T_i}\n")


def plot_cogging_torque():
    dir_data = Path(__file__).parent / "data"
    dir_media = Path(__file__).parent / "media"

    plt.figure()

    theta = []
    T = []
    with open(dir_data / "cogging_torque.csv") as f:
        for line in f:
            line = [float(li) for li in line.strip().split(",")]
            T.append(line.pop() * 8)
            theta.append(line.pop())

    plt.plot(theta, T, "b-")
    plt.grid()
    plt.xlabel("Rotor angle [°]")
    plt.ylabel("Torque [Nm]")
    plt.savefig(dir_media / "cogging_torque.png", dpi=550, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    plt.style.use(["default", "fast"])
    # calculate_torque(I0=50)
    # calculate_torque(I0=100)
    # calculate_torque(I0=150)
    # calculate_torque(I0=200)
    # calculate_torque(I0=250)

    calculate_cogging_torque()
    plot_cogging_torque()