from itertools import product

import matplotlib.pyplot as plt
import numpy as np

from digital_twin_distiller import ModelDir
from numpy import linspace
import pandas as pd
import json

ModelDir.set_base(__file__)

range_a0 = 0.5
range_a1 = 3.5
nsteps_a = 31

range_b0 = 0.0
range_b1 = 3.0
nsteps_b = 31

range_c0 = 0.0
range_c1 = 45
nsteps_c = 25

range_a = linspace(range_a0, range_a1, nsteps_a)
range_b = linspace(range_b0, range_b1, nsteps_b)
range_c = linspace(range_c0, range_c1, nsteps_c)

prod = list(product(range_a, range_b, range_c))
range_prod = linspace (0, len(prod), len(prod)+1)
prod1 = list(product(range_a, range_b))

f = open(ModelDir.DATA / f'locked_rotorv2.json')
torque = json.load(f)

res = {}
res["earheight"] = [(prod[i])[0] for i in range(len(prod))]
res["aslheight"] = [(prod[i])[1] for i in range(len(prod))]
res["rotorangle"] = [(prod[i])[2] for i in range(len(prod))]
res["torque"] = [(torque["Torque"])[i] for i in range(len(prod))]
res = pd.DataFrame(res)

switch = 1
if switch == 0:
    t = [[] for i in range(len(prod1))]
    a = 0
    b = 0
    while a < len(prod1):
        t[a] = [(torque["Torque"])[i] for i in range(b+0, b+25)]
        a = a + 1
        b = b + 25

    tmax = [[] for i in range(len(prod1))]
    tmid = [[] for i in range(len(prod1))]
    inmid = [[] for i in range(len(prod1))]
    tpeak = [[0 for i in range(nsteps_c)] for i in range(len(prod1))]
    inpeak = [[0 for i in range(nsteps_c)] for i in range(len(prod1))]
    tmin = [[] for i in range(len(prod1))]

    for i in range(len(prod1)):
        tmax[i] = max(t[i])
        tmin[i] = min(t[i])
    tabsmax = max(abs(x) for x in tmax)
    tabsmin = min(abs(x) for x in tmin)

    for a in range(len(prod1)):
        for b in range(1, nsteps_c-1):
            if (t[a])[b] > (t[a])[b-1] and (t[a])[b] > (t[a])[b+1]:
                (tpeak[a])[b] = round((t[a])[b], 5)
            elif (t[a])[b] < (t[a])[b-1] and (t[a])[b] < (t[a])[b+1]:
                (tpeak[a])[b] = round((t[a])[b], 5)
            else:
                pass
            if (tpeak[a])[b] == 0:
                pass
            else:
                (inpeak[a])[b] = round(b / 25 * 45, 1)
        inpeak[a] = [i for i in inpeak[a] if i != 0]
        tpeak[a] = [i for i in tpeak[a] if i != 0]

    case = {"earheight": [(prod1[i])[0] for i in range(len(prod1))],
            "aslheight": [(prod1[i])[1] for i in range(len(prod1))],
            "coggingtorque": t,
            "torquepeak": tpeak,
            "peakindex": inpeak}
    case = pd.DataFrame(case)

elif switch == 1:
    t = [[] for i in range(nsteps_a)]
    a = 0
    b = 0
    while a < nsteps_a:
        t[a] = [(torque["Torque"])[i] for i in range(b + 0, b + 25)]
        a = a + 1
        b = b + 650
    print(t)

    tmax = [[] for i in range(nsteps_a)]
    tmid = [[] for i in range(nsteps_a)]
    inmid = [[] for i in range(nsteps_a)]
    tpeak = [[0 for i in range(nsteps_c)] for i in range(nsteps_a)]
    inpeak = [[0 for i in range(nsteps_c)] for i in range(nsteps_a)]
    tmin = [[] for i in range(nsteps_a)]

    for i in range(nsteps_a):
        tmax[i] = max(t[i])
        tmin[i] = min(t[i])
    tabsmax = max(abs(x) for x in tmax)
    tabsmin = min(abs(x) for x in tmin)

    for a in range(nsteps_a):
        for b in range(1, nsteps_c - 1):
            if (t[a])[b] > (t[a])[b - 1] and (t[a])[b] > (t[a])[b + 1]:
                (tpeak[a])[b] = round((t[a])[b], 5)
            elif (t[a])[b] < (t[a])[b - 1] and (t[a])[b] < (t[a])[b + 1]:
                (tpeak[a])[b] = round((t[a])[b], 5)
            else:
                pass
            if (tpeak[a])[b] == 0:
                pass
            else:
                (inpeak[a])[b] = round(b / 25 * 45, 1)
        inpeak[a] = [i for i in inpeak[a] if i != 0]
        tpeak[a] = [i for i in tpeak[a] if i != 0]

    case = {"earheight": range_a,
            "coggingtorque": t,
            "torquepeak": tpeak,
            "peakindex": inpeak}
    case = pd.DataFrame(case)

elif switch == 3:
    t = [[] for i in range(nsteps_a)]
    a = 0
    b = 0
    while a < nsteps_a:
        t[a] = [(torque["Torque"])[i] for i in range(b + 0, b + 76)]
        t[a] = [round((t[a])[i], 3) for i in range(len(range_c))]
        a = a + 1
        b = b + 2356

    tmax = [[] for i in range(nsteps_a)]
    tmin = [[] for i in range(nsteps_a)]
    inmax = [[] for i in range(nsteps_a)]
    inmin = [[] for i in range(nsteps_a)]

    for i in range(nsteps_a):
        tmax[i] = max(t[i])
        tmin[i] = min(t[i])
        inmax[i] = np.multiply(t[i].index(tmax[i]), 0.1)
        inmin[i] = np.multiply(t[i].index(tmin[i]), 0.1)
        if inmin[i] == 7.5:
            inmin[i] = None
            tmin[i] = None
        else:
            inmin[i] = inmin[i]
            tmin[i] = tmin[i]
    case = {"earheight": range_a,
            "coggingtorque": t,
            "inmaxpeak": inmax,
            "inminpeak": inmin,
            "tmaxpeak": tmax,
            "tminpeak": tmin}
    case = pd.DataFrame(case)

    case.to_pickle(ModelDir.DATA / "df_cogging.pkl")

elif switch == 4:
    t = [[] for i in range(len(prod1))]
    a = 0
    b = 0
    while a < len(prod1):
        t[a] = [(torque["Torque"])[i] for i in range(b + 0, b + 76)]
        t[a] = [round((t[a])[i], 3) for i in range(len(range_c))]
        a = a + 1
        b = b + 76

    tmax = [[] for i in range(len(prod1))]
    tmin = [[] for i in range(len(prod1))]
    inmax = [[] for i in range(len(prod1))]
    inmin = [[] for i in range(len(prod1))]
    tdelta1 = [[] for i in range(len(prod1))]
    tdelta2 = [[] for i in range(len(prod1))]

    for i in range(len(prod1)):
        tmax[i] = max(t[i])
        tmin[i] = min(t[i])
        inmax[i] = np.multiply(t[i].index(tmax[i]), 0.1)
        inmin[i] = np.multiply(t[i].index(tmin[i]), 0.1)
        tdelta1[i] = abs(tmax[i] - tmin[i])
        tdelta2[i] = tmax[i] - abs(tmin[i])

    case = {"earheight": [(prod1[i])[0] for i in range(len(prod1))],
            "aslheight": [(prod1[i])[1] for i in range(len(prod1))],
            "coggingtorque": t,
            "inmaxpeak": inmax,
            "inminpeak": inmin,
            "tmaxpeak": tmax,
            "tminpeak": tmin,
            "tdelta1": tdelta1,
            "tdelta2": tdelta2
            }
    case = pd.DataFrame(case)

    x = case[case["tdelta2"] == min(case["tdelta2"], key=abs)]
    print(x)

    case.to_pickle(ModelDir.DATA / "df_cogging.pkl")

else:
    pass

switch = 3
if switch == 0:
    fig, ax1 = plt.subplots()
    plt.title("Peak of the cogging torque and the rotor angle")
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
    ax1.set_ylabel('Peak cogging torque [Nm]', c="b")
    ax2.set_ylabel('Rotor angle [°]', c="r")
    plt.show()

    a = 0
    b = 1
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
    plt.show()

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

    c = 0
    a = 0
    b = 961
    for i in range(a, b):
        ax = plt.axes(projection='3d')
        zdata = (case["torquepeak"].loc[range(a, b)]).tolist()
        xdata = (case["earheight"].loc[range(a, b)]).tolist()
        ydata = (case["aslheight"].loc[range(a, b)]).tolist()
        ax.scatter3D(xdata, ydata, zdata)
    plt.show()

elif switch == 3:

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
    b = 31
    for i in range(a, b, 5):
        sub1.plot(range(33, 40), [((case['coggingtorque'])[i])[a] for a in range(15, 22)])
        sub2.plot(range_c, (case["coggingtorque"])[i], label=("A=" + str(0.5+0*0.1) + "mm"  + "," + " C=" + str(i*0.1) + "mm"))
    sub1.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    sub1.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    sub1.minorticks_on()
    sub2.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    sub2.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    sub2.minorticks_on()
    sub1.set_xticklabels([0, 28, 30, 32, 34, 36, 38, 40])
    plt.xlabel("Electrical angle [deg]", fontsize=10)
    plt.ylabel("Torque [Nm]", fontsize=10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(bbox_to_anchor=(0.64, 1.225), fontsize=8.5)
    #plt.savefig(ModelDir.MEDIA / "cogging.png", bbox_inches="tight", dpi=650)
    plt.show()


else:
    pass