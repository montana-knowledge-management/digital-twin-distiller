from itertools import product

import matplotlib.pyplot as plt
import numpy as np

from digital_twin_distiller import ModelDir
from numpy import linspace
import pandas as pd
from matplotlib.lines import Line2D

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
range_prod = linspace(0, len(prod), len(prod)+1)
prod1 = list(product(range_a, range_b))

range_c = linspace(range_c0, range_c1, nsteps_c)

case = pd.read_pickle(ModelDir.DATA / "df_cogging.pkl")

switch = 0
if switch == 0:
    c = 16
    a = 0
    b = 806
    fig, ax1 = plt.subplots(figsize=(6, 4))
    #ax2 = ax1.twinx()
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["tmaxpeak"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="#B90276")
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["tminpeak"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c='#005691')
    #for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tdelta3"].loc[range(a, b)]):
        #ax2.scatter(xe, ye, c="b")
    #for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tdelta3"].loc[range(a, b)]):
       #ax2.scatter(xe, ye, c="g")
    ax1.set_xlabel('"A" paraméter [mm]', fontsize=12)
    ax1.set_ylabel('Nyomaték [Nm]', fontsize=12)
    #ax2.set_ylabel('Torque [Nm]', fontsize=10, c='r')
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    #ax2.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color="#B90276", label= "maximum"),
              Line2D([0], [0], marker="o", color='#005691', label= "minimum")]
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(handles=legend, fontsize = 10)
    plt.savefig(ModelDir.MEDIA / "mait2.png", bbox_inches="tight", dpi=650)
    plt.show()

    fig, ax1 = plt.subplots(figsize=(6, 4))
    # ax2 = ax1.twinx()
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["inmaxpeak"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="#008ECF")
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["inminpeak"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c='#78BE20')
    # for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tdelta3"].loc[range(a, b)]):
    # ax2.scatter(xe, ye, c="b")
    # for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tdelta3"].loc[range(a, b)]):
    # ax2.scatter(xe, ye, c="g")
    ax1.set_xlabel('"A" paraméter" [mm]', fontsize=12)
    ax1.set_ylabel('Terhelési szög [°]', fontsize=12)
    # ax2.set_ylabel('Torque [Nm]', fontsize=10, c='r')
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    # ax2.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color="#008ECF", label= "maximum"),
              Line2D([0], [0], marker="o", color='#78BE20', label= "minimum")]
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(handles=legend, fontsize = 10)
    plt.savefig(ModelDir.MEDIA / "mait3.png", bbox_inches="tight", dpi=650)
    plt.show()

    x = [[] for i in range(805)]
    rev = case.copy()
    rev = rev.loc[::-1].reset_index(drop=True)
    rev["coggingtorque"] = list(reversed(rev["coggingtorque"]))
    for i in range(805):
        x[i] = [a * -1 for a in (rev["coggingtorque"])[i]]
    rev = {'x': x}
    rev2 = case.copy()
    for i in range(805):
        x[i] = [a * -1 for a in (rev2["coggingtorque"])[i]]
    rev2 = {'x': x}
    c = 20
    a = c*31
    b = c*31+1
    fig = plt.figure(figsize=(6, 4))
    for i in range(a, b):
        plt.plot(range_c * 4, (case["coggingtorque"])[i], c='b')
        plt.plot((7.5-range_c) * 4, (rev["x"])[i], c='b')
        plt.plot((15.05-range_c) * 4, (rev2["x"])[i], c='b')
        plt.plot((7.55+range_c) * 4, (case["coggingtorque"])[i], c='b')
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Electrical angle [deg]", fontsize=10)
    plt.ylabel("Torque [Nm]", fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    #plt.savefig(ModelDir.MEDIA / "cogging2.png", bbox_inches="tight", dpi=650)
    #plt.show()

    a = 620
    b = 651
    for i in range(a, b):
        plt.plot(range_c, (case["coggingtorque"])[i])
        plt.plot((case["inmaxpeak"])[i], (case["tmaxpeak"])[i], "x")
        plt.plot((case["inminpeak"])[i], (case["tminpeak"])[i], "x")
    #plt.show()

elif switch == 1:
    c = 19
    a = 31 * c
    b = 31 * (c+1)
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
    plt.show()

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
        ydata = (case["inmaxpeak"].loc[range(a, b)]).tolist()
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

elif switch == 4:
    a = 0
    b = 806
    c = 5

    c = 31 * c
    t = [[]*26]
    colors = ["#B90276", '#50237F', '#005691', "#008ECF", '#78BE20', "#006249", '#525F6B', '#00A8B0', '#000']
    fig = plt.subplots(figsize=(6, 4))
    for i, j in zip(range(a, b, c), range(9)):
        t = (case["coggingtorque"])[i]
        t = np.multiply(t,-1)
        tt = t[0]
        x = [t[d] - 2*tt for d in range(76)]


        plt.plot(range_c*4, (case["coggingtorque"])[i], color = colors[j], label="A=" + str(0.5+i/310) + "mm")
        plt.plot(7.5*4-range_c*4, x, color = colors[j])
        plt.xlabel('Terhelési szög [°]', fontsize=12)
        plt.ylabel('Nyomaték [Nm]', fontsize=12)
        plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
        plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
        plt.minorticks_on()
        plt.xticks(np.arange(0, 32, step=3))
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=10)
    plt.savefig(ModelDir.MEDIA / "mait1.png", bbox_inches="tight", dpi=650)
    plt.show()