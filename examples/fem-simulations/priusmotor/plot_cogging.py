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

range_c0 = 3.75
range_c1 = 7.5
nsteps_c = 76

range_a = linspace(range_a0, range_a1, nsteps_a)
range_b = linspace(range_b0, range_b1, nsteps_b)
range_c = linspace(range_c0, range_c1, nsteps_c)

prod = list(product(range_a, range_b, range_c))
range_prod = linspace (0, len(prod), len(prod)+1)
prod1 = list(product(range_a, range_b))

range_c = linspace(range_c0, range_c1, nsteps_c)

case = pd.read_pickle(ModelDir.DATA / "df_cogging.pkl")

switch = 0
if switch == 0:
    fig, ax1 = plt.subplots()
    plt.title("Peak of the cogging torque and the rotor angle")
    ax2 = ax1.twinx()
    for xe, ye in zip(case["earheight"], case["tdelta1"]):
        ax1.scatter(xe, ye, c="blue")
    #for xe, ye in zip(case["earheight"], case["tminpeak"]):
        #ax1.scatter(xe, ye, c="blue")
    for xe, ye in zip(case["earheight"], case["tdelta2"]):
        ax2.scatter(xe, ye, c="g")
    #for xe, ye in zip(case["earheight"], case["inminpeak"]):
       #ax2.scatter(xe, ye, c="red")
    ax1.set_xlabel('Variable parameter [mm]')
    ax1.set_ylabel('Peak cogging torque [Nm]', c="b")
    ax2.set_ylabel('Rotor angle [°]', c="r")
    plt.show()

    a = 589
    b = 620
    for i in range(a, b):
        plt.plot(range_c, (case["coggingtorque"])[i], lw=2)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("rotorangle [deg]")
    plt.ylabel("Torque [Nm]")
    plt.show()

    a = 0
    b = 31
    for i in range(a, b):
        plt.plot(range_c, (case["coggingtorque"])[i])
        plt.plot((case["inmaxpeak"])[i], (case["tmaxpeak"])[i], "x")
        plt.plot((case["inminpeak"])[i], (case["tminpeak"])[i], "x")
    #plt.show()

elif switch == 1:
    c = 19
    a = 0
    b = 805
    fig, ax1 = plt.subplots()
    plt.title("Peak of the cogging torque and the rotor angle")
    ax2 = ax1.twinx()
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tdelta1"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="blue")
    #for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tminpeak"].loc[range(a, b)]):
        #ax1.scatter(xe, ye, c="red")
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tdelta2"].loc[range(a, b)]):
        ax2.scatter(xe, ye, c="green")
    #for xe, ye in zip(case["tdelta1"].loc[range(a, b)], case["inminpeak"].loc[range(a, b)]):
        #ax2.scatter(xe, ye, c="red")
    ax1.set_xlabel('Torque difference 2 [Nm]')
    ax1.set_ylabel('Torque difference 1 [Nm]', c="b")
    #ax2.set_ylabel('Torque difference 2 [Nm]', c="r")
    #plt.show()

    a = c * 31
    b = (c + 1) * 31
    for i in range(a, b):
        plt.plot(range_c, (case["coggingtorque"])[i])
        plt.plot((case["inmaxpeak"])[i], (case["tmaxpeak"])[i], "x")
        plt.plot((case["inminpeak"])[i], (case["tminpeak"])[i], "x")
    #plt.show()

elif switch == 2:

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection='3d')
    c = 0
    a = 0
    b = 806
    for i in range(a, b):
        zdata = (case["tmaxpeak"].loc[range(a, b)]).tolist()
        xdata = (case["earheight"].loc[range(a, b)]).tolist()
        ydata = (case["aslheight"].loc[range(a, b)]).tolist()
        ax.scatter3D(xdata, ydata, zdata)
        ax.set_xlabel('Parameter A [mm]', fontsize=10)
        ax.set_ylabel('Parameter C [mm]', fontsize=10)
        ax.set_zlabel('Torque [Nm]', fontsize=10)
        ax.minorticks_on()
        ax.tick_params(labelsize=10)
    plt.savefig(ModelDir.MEDIA / "cogging3d.png", bbox_inches="tight", dpi=650)
    plt.show()

    case['minmax'] = list(zip(case["tminpeak"], case["tmaxpeak"]))

    #fig = plt.figure(figsize=(6, 6))
    #ax = fig.add_subplot(projection='3d')
    #c = 0
    #a = 0
    #b = 806
    #for i in range(b):
        #zdata = ((case["minmax"])[i])
        #xdata = ((case["earheight"])[i])
        #ydata = ((case["aslheight"])[i])
        #ax.scatter3D(xdata, ydata, zdata)
        #ax.set_xlabel('Parameter A [mm]', fontsize=10)
        #ax.set_ylabel('Parameter C [mm]', fontsize=10)
        #ax.set_zlabel('Torque [Nm]', fontsize=10)
        #ax.minorticks_on()
        #ax.tick_params(labelsize=10)
    #plt.show()