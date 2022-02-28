from itertools import product

import matplotlib.pyplot as plt
from digital_twin_distiller import ModelDir
from numpy import linspace
import pandas as pd

ModelDir.set_base(__file__)

range_a0 = 0.5
range_a1 = 3.0
nsteps_a = 26

range_b0 = 0.0
range_b1 = 3.0
nsteps_b = 31

range_c0 = 0.0
range_c1 = 45
nsteps_c = 91

range_a = linspace(range_a0, range_a1, nsteps_a)
range_b = linspace(range_b0, range_b1, nsteps_b)
range_c = linspace(range_c0, range_c1, nsteps_c)

prod = list(product(range_a, range_b, range_c))
range_prod = linspace (0, len(prod), len(prod)+1)
prod1 = list(product(range_a, range_b))

range_c = linspace(range_c0, range_c1, nsteps_c)

case = pd.read_pickle(ModelDir.DATA / "df_locked.pkl")

switch = 2
if switch == 0:
    fig, ax1 = plt.subplots()
    plt.title("Peak of the torque and the rotor angle")
    ax2 = ax1.twinx()
    for xe, ye in zip(case["earheight"], case["tmaxpeak"]):
        ax1.scatter(xe, ye, c="blue")
    for xe, ye in zip(case["earheight"], case["tminpeak"]):
        ax1.scatter(xe, ye, c="blue")
    for xe, ye in zip(case["earheight"], case["inmaxpeak"]):
        ax2.scatter(xe, ye, c="red")
    for xe, ye in zip(case["earheight"], case["inminpeak"]):
        ax2.scatter(xe, ye, c="red")
    ax1.set_xlabel('Variable parameter [mm]')
    ax1.set_ylabel('Peak torque [Nm]', c="b")
    ax2.set_ylabel('Rotor angle [°]', c="r")
    plt.show()

    a = 0
    b = 1
    for i in range(a, b):
        plt.plot(range_c, (case["torque"])[i], lw=2)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("rotorangle [deg]")
    plt.ylabel("torque [Nm]")
    plt.show()

    a = 0
    b = 31
    for i in range(a, b):
        plt.plot(range_c, (case["torque"])[i])
        plt.plot((case["inmaxpeak"])[i], (case["tmaxpeak"])[i], "x")
        plt.plot((case["inminpeak"])[i], (case["tminpeak"])[i], "x")
    plt.show()

elif switch == 1:
    a = 0
    b = 1
    fig, ax1 = plt.subplots()
    plt.title("Peak of the torque and the rotor angle")
    ax2 = ax1.twinx()
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["torquepeak"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="blue")
    #for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tminpeak"].loc[range(a, b)]):
        #ax1.scatter(xe, ye, c="blue")
    #for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["inmaxpeak"].loc[range(a, b)]):
        #ax2.scatter(xe, ye, c="red")
    #for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["inminpeak"].loc[range(a, b)]):
        #ax2.scatter(xe, ye, c="red")
    ax1.set_xlabel('Variable parameter [mm]')
    ax1.set_ylabel('Peak torque [Nm]', c="b")
    ax2.set_ylabel('Rotor angle [°]', c="r")
    plt.show()

elif switch == 2:
    a = 0
    b = 1
    for i in range(a, b):
        plt.plot(list(range(nsteps_c)), (case["torque"])[i])
        plt.plot((case["maxpeak"])[i], [((case["torque"])[i])[a] for a in (case["maxpeak"])[i]], "x")
        plt.plot((case["minpeak"])[i], [((case["torque"])[i])[a] for a in (case["minpeak"])[i]], "x")
    plt.show()