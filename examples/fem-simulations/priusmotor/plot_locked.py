from itertools import product
import numpy as np
import matplotlib.pyplot as plt
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

switch = -1
if switch == -2:
    fig = plt.figure(figsize=(6, 4))
    a = 0
    b = 31 * 26
    range_c = np.multiply(range_c, 4)
    colors = ["#B90276", '#50237F', '#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    for i,j, e in zip(range(a, b, 31*5), range(5, 31, 5), range(9)):
        plt.plot(range(110, 158, 2), [((case['torque'])[i])[a] for a in range(55, 79)], label=("A=" + str(j*0.1) + "mm"), color=colors[e])
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(110, 160, step=4))
    plt.xlabel("Terhelési szög [°]", fontsize=12)
    plt.ylabel("Nyomaték [Nm]", fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=10)
    plt.savefig(ModelDir.MEDIA / "mait4", bbox_inches="tight", dpi=650)
    plt.show()

if switch == -1:
    fig = plt.figure(figsize=(6, 4))
    a = 31 * 15
    b = 31 * 16
    range_c = np.multiply(range_c, 4)
    colors = ["#B90276", '#50237F', '#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    for i,j, e in zip(range(a, b, 6), range(0, 31, 6), range(9)):
        plt.plot(range(110, 158, 2), [((case['torque'])[i])[a] for a in range(55, 79)], label=("A=" + str(2.1) + "mm" + "," + " C=" + str(round(j*0.1, 2)) + "mm"), color=colors[e])
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Terhelési szög[°]", fontsize=12)
    plt.ylabel("Nyomaték [Nm]", fontsize=12)
    plt.xticks(np.arange(110, 160, step=4))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=10)
    plt.savefig(ModelDir.MEDIA / "mait5.png", bbox_inches="tight", dpi=650)
    plt.show()

if switch == 0:
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection='3d')
    c = 0
    a = 0
    b = 806
    for i in range(a, b):
        zdata = (case["i4"].loc[range(a, b)]).tolist()
        xdata = (case["earheight"].loc[range(a, b)]).tolist()
        ydata = (case["aslheight"].loc[range(a, b)]).tolist()
        ax.scatter3D(xdata, ydata, zdata)
        ax.set_xlabel('Parameter A [mm]', fontsize=10)
        ax.set_ylabel('Parameter C [mm]', fontsize=10)
        ax.set_zlabel('Electric angle [deg]', fontsize=10)
        ax.minorticks_on()
        ax.view_init(elev=20, azim=280)
        ax.set_ylim(3, 0)
        ax.set_ylim(3, 0)
        ax.tick_params(labelsize=10)
    plt.savefig(ModelDir.MEDIA / "i4_locked3d.png", bbox_inches="tight", dpi=650)
    #plt.show()

if switch == 1:
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection='3d')
    c = 0
    a = 0
    b = 806
    for i in range(a, b):
        zdata = (case["i4"].loc[range(a, b)]).tolist()
        xdata = (case["earheight"].loc[range(a, b)]).tolist()
        ydata = (case["aslheight"].loc[range(a, b)]).tolist()
        ax.scatter3D(xdata, ydata, zdata)
        ax.set_xlabel('Parameter A [mm]', fontsize=10)
        ax.set_ylabel('Parameter C [mm]', fontsize=10)
        ax.set_zlabel('Electric angle [deg]', fontsize=10)
        ax.minorticks_on()
        ax.view_init(elev=20, azim=280)
        ax.set_ylim(3, 0)
        ax.set_ylim(3, 0)
        ax.tick_params(labelsize=10)
    plt.savefig(ModelDir.MEDIA / "i4_locked3d.png", bbox_inches="tight", dpi=650)
    #plt.show()

elif switch == 2:
    fig = plt.figure(figsize=(6, 4))
    plt.subplots_adjust(bottom=0.1, left=0.1, top=0.9, right=0.9)
    sub1 = fig.add_subplot(2, 3, (1,2))
    sub2 = fig.add_subplot(2, 3,(4,6))
    #a = 14
    #z = np.polyfit(range_c, case["coggingtorque"][0], a)
    #predict = np.poly1d(z)
    #y = predict(range_c)
    #plt.plot(range_c, y, c="r", linestyle='--')
    a = 31 * 15
    b = 31 * 16
    range_c = np.multiply(range_c, 4)
    for i, j in zip(range(a, b, 5), range(0, 31, 5)):
        sub1.plot(range(110, 162, 2), [((case['torque'])[i])[a] for a in range(55, 81)])
        sub2.plot(range_c, (case["torque"])[i], label=("A=" + str(2.1) + "mm"  + "," + " C=" + str(j*0.1) + "mm"))
    sub1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    sub1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    sub1.minorticks_on()
    sub2.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    sub2.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    sub2.minorticks_on()
    plt.xlabel("Electrical angle [deg]", fontsize=10)
    plt.ylabel("Torque [Nm]", fontsize=10)
    #sub1.set_xticks(np.arange(110, 161, step=10), fontsize=10)
    sub2.set_xticks(np.arange(0, 200, step=20), fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(bbox_to_anchor=(0.64, 1.225), fontsize=8.5)
    plt.savefig(ModelDir.MEDIA / "PEMC_locked.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 3:
    a = 0
    b = 26
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["tmaxpeak"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="r")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10, c="b")
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='r', label="tmax")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "tmax1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 4:
    a = 0
    b = 26
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["tminpeak"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "tmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 5:
    a = 0
    b = 26
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["inmaxpeak"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Electrical angle [deg]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="inmax")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "inmax1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 6:
    a = 0
    b = 26
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["inminpeak"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 7:
    a = 0
    b = 26
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["t2"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    #plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 8:
    a = 0
    b = 26
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["t3"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    #plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 9:
    a = 0
    b = 26
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["i2"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    # plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 10:
    a = 0
    b = 26
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["earheight"].loc[range(a, b)], case["i3"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    # plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

if switch == 11:
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
        #plt.gca().invert_yaxis()
        #ax.view_init(elev=20., azim=140)
    #plt.savefig(ModelDir.MEDIA / "cogging3d.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 12:
    fig = plt.figure(figsize=(6, 4))
    plt.subplots_adjust(bottom=0.1, left=0.1, top=0.9, right=0.9)
    sub1 = fig.add_subplot(2, 3, (1,2))
    sub2 = fig.add_subplot(2, 3,(4,6))
    #a = 14
    #z = np.polyfit(range_c, case["coggingtorque"][0], a)
    #predict = np.poly1d(z)
    #y = predict(range_c)
    #plt.plot(range_c, y, c="r", linestyle='--')
    a = 0
    b = 26
    for i in range(a, b, 5):
        sub1.plot(range(55, 80), [((case['torque'])[i])[a] for a in range(55, 80)])
        sub2.plot(range_c, (case["torque"])[i], label=("A=" + str(0.5+0*0.1) + "mm"  + "," + " C=" + str(i*0.1) + "mm"))
    sub1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    sub1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    sub1.minorticks_on()
    sub2.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    sub2.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    sub2.minorticks_on()
    plt.xlabel("Electrical angle [deg]", fontsize=10)
    plt.ylabel("Torque [Nm]", fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(bbox_to_anchor=(0.64, 1.225), fontsize=8.5)
    plt.savefig(ModelDir.MEDIA / "cogging.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 13:
    a = 0
    b = 31
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tmax"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="r")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10, c="b")
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='r', label="tmax")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "tmax1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 14:
    a = 0
    b = 31
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tmin"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "tmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 15:
    a = 0
    b = 31
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["inmax"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Electrical angle [deg]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="inmax")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "inmax1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 16:
    a = 0
    b = 31
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["inmin"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 17:
    a = 0
    b = 31
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["t2"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    #plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 18:
    a = 0
    b = 31
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["t3"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    #plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 19:
    a = 0
    b = 31
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["i2"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    # plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 20:
    a = 0
    b = 31
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["i3"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    # plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

if switch == 21:
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
        #plt.gca().invert_yaxis()
        #ax.view_init(elev=20., azim=140)
    #plt.savefig(ModelDir.MEDIA / "cogging3d.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 22:
    fig = plt.figure(figsize=(6, 4))
    plt.subplots_adjust(bottom=0.1, left=0.1, top=0.9, right=0.9)
    sub1 = fig.add_subplot(2, 3, (1,2))
    sub2 = fig.add_subplot(2, 3,(4,6))
    #a = 14
    #z = np.polyfit(range_c, case["coggingtorque"][0], a)
    #predict = np.poly1d(z)
    #y = predict(range_c)
    #plt.plot(range_c, y, c="r", linestyle='--')
    a = 0
    b = 26
    for i in range(a, b, 5):
        sub1.plot(range(55, 80), [((case['torque'])[i])[a] for a in range(55, 80)])
        sub2.plot(range_c, (case["torque"])[i], label=("A=" + str(0.5+0*0.1) + "mm"  + "," + " C=" + str(i*0.1) + "mm"))
    sub1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    sub1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    sub1.minorticks_on()
    sub2.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    sub2.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    sub2.minorticks_on()
    plt.xlabel("Electrical angle [deg]", fontsize=10)
    plt.ylabel("Torque [Nm]", fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(bbox_to_anchor=(0.64, 1.225), fontsize=8.5)
    plt.savefig(ModelDir.MEDIA / "cogging.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 23:
    a = 0
    b = 806
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tmax"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="r")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10, c="b")
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='r', label="tmax")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "tmax1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 24:
    a = 0
    b = 806
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tmin"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "tmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 25:
    a = 0
    b = 806
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["inmax"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Electrical angle [deg]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="inmax")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "inmax1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 26:
    a = 0
    b = 806
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["inmin"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 27:
    a = 0
    b = 806
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["t2"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    #plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 28:
    a = 0
    b = 806
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["t3"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    #plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 29:
    a = 0
    b = 806
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["i2"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    # plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 30:
    a = 0
    b = 806
    fig, ax1 = plt.subplots(figsize=(6, 4))
    for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["i3"].loc[range(a, b)]):
        ax1.scatter(xe, ye, c="b")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10)
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='b', label="tmin")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    # plt.savefig(ModelDir.MEDIA / "inmin1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 31:
    tt = [[0]*21 for i in range(806)]
    a = 31 * 16
    b = 31 * 17
    for i in range(a, b):
        tt[i] = ((case["torque"])[i])[5:10]
        plt.plot(range(5, 10), tt[i])
        maxi = max(tt[i])
        print(maxi)
    plt.show()
elif switch == 99:
    a = 0
    b = 26
    fig, ax1 = plt.subplots(figsize=(6, 4))
    #ax2 = ax1.twinx()
    for xe, ye in zip(case["earheight"], (case["tmaxpeak"])[0]):
        ax1.scatter(xe, ye, c="r")
    #for xe, ye in zip(case["earheight"].loc[range(a, b)], case["tminpeak"].loc[range(a, b)]):
        #ax1.scatter(xe, ye, c="b")
    #for xe, ye in zip(case["earheight"].loc[range(a, b)], case["tminpeak"].loc[range(a, b)]):
        #ax2.scatter(xe, ye, c="b")
    # for xe, ye in zip(case["aslheight"].loc[range(a, b)], case["tdelta3"].loc[range(a, b)]):
    # ax2.scatter(xe, ye, c="g")
    ax1.set_xlabel('Parameter A [mm]', fontsize=10)
    ax1.set_ylabel('Torque [Nm]', fontsize=10, c="b")
    #ax2.set_ylabel('Torque [Nm]', fontsize=10, c='r')
    ax1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax1.minorticks_on()
    #ax2.minorticks_on()
    legend = [Line2D([0], [0], marker="o", color='r', label="tmax"),
              Line2D([0], [0], marker="o", color='b', label=u"\u0394" + "t2")]
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(handles=legend)
    #plt.savefig(ModelDir.MEDIA / "delta.png", bbox_inches="tight", dpi=650)
    plt.show()