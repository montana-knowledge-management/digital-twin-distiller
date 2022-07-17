import matplotlib.pyplot as plt
import numpy as np
from numpy import linspace
from itertools import product
import json
import pandas as pd
from digital_twin_distiller import ModelDir
from scipy.signal import find_peaks

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
range_prod = linspace(0, len(prod), len(prod)+1)
prod1 = list(product(range_a, range_b))
prod2 = list(product(range_a, range_c))
prod3 = list(product(range_b, range_c))

f = open(ModelDir.DATA / f'locked_rotor.json')
torque = json.load(f)

switch = 2
if switch == 1:
    t = [[] for i in range(len(range_a))]
    a = 0
    b = 0
    while a < len(range_a):
        t[a] = [(torque["Torque"])[i] for i in range(b + 0, b + 91)]
        t[a] = [round((t[a])[i], 3) for i in range(len(range_c))]
        a = a + 1
        b = b + 2821

    tmax = [[] for i in range(len(range_a))]
    tmin = [[] for i in range(len(range_a))]
    inmax = [[] for i in range(len(range_a))]
    inmin = [[] for i in range(len(range_a))]
    tempmax = [[] for i in range(len(range_a))]
    tempmin = [[] for i in range(len(range_a))]
    tmaxpeaks = [[0] * 2 for i in range(len(range_a))]
    inmaxpeaks = [[0] * 2 for i in range(len(range_a))]
    tminpeaks = [[0] * 2 for i in range(len(range_a))]
    inminpeaks = [[0] * 2 for i in range(len(range_a))]
    t1 = [[] for i in range(len(range_a))]
    t2 = [[] for i in range(len(range_a))]
    t3 = [[] for i in range(len(range_a))]
    t4 = [[] for i in range(len(range_a))]
    i1 = [[] for i in range(len(range_a))]
    i2 = [[] for i in range(len(range_a))]
    i3 = [[] for i in range(len(range_a))]
    i4 = [[] for i in range(len(range_a))]

    for i in range(len(range_a)):
        tmax[i] = max(t[i])
        tmin[i] = min(t[i])
        inmax[i] = np.multiply(t[i].index(tmax[i]), 2)
        inmin[i] = np.multiply(t[i].index(tmin[i]), 2)
        tempmax[i], _ = find_peaks(t[i])
        tempmin[i], _ = find_peaks(np.multiply(t[i], -1))
        inmaxpeaks[i] = np.multiply(tempmax[i], 2)
        inminpeaks[i] = np.multiply(tempmin[i], 2)
        for j in range(len(tempmax[i])):
            (tmaxpeaks[i])[j] = (t[i])[(tempmax[i])[j]]
            (tminpeaks[i])[j] = (t[i])[(tempmin[i])[j]]
    for i in range(len(range_a)):
        t1[i] = (tminpeaks[i])[0]
        t2[i] = (tmaxpeaks[i])[0]
        t3[i] = (tminpeaks[i])[1]
        t4[i] = (tmaxpeaks[i])[1]
        i1[i] = (inminpeaks[i])[0]
        i2[i] = (inmaxpeaks[i])[0]
        i3[i] = (inminpeaks[i])[1]
        i4[i] = (inmaxpeaks[i])[1]

    case = {"earheight": range_a,
            "torque": t,
            "inmax": inmax,
            "inmin": inmin,
            "tmax": tmax,
            "tmin": tmin,
            "inmaxpeaks": inmaxpeaks,
            "inminpeaks": inminpeaks,
            "tmaxpeaks": tmaxpeaks,
            "tminpeaks": tminpeaks,
            "t1": t1,
            "t2": t2,
            "t3": t3,
            "t4": t4,
            "i1": i1,
            "i2": i2,
            "i3": i3,
            "i4": i4}
    case = pd.DataFrame(case)
    case.to_pickle(ModelDir.DATA / "df_locked.pkl")

elif switch == 2:
    t = [[] for i in range(len(prod1))]
    a = 0
    b = 0
    while a < len(prod1):
        t[a] = [(torque["Torque"])[i] for i in range(b + 0, b + 91)]
        t[a] = [round((t[a])[i], 3) for i in range(len(range_c))]
        a = a + 1
        b = b + 91

    tmax = [[] for i in range(len(prod1))]
    tmin = [[] for i in range(len(prod1))]
    inmax = [[] for i in range(len(prod1))]
    inmin = [[] for i in range(len(prod1))]
    tempmax = [[] for i in range(len(prod1))]
    tempmin = [[] for i in range(len(prod1))]
    tmaxpeaks = [[0] * 2 for i in range(len(prod1))]
    inmaxpeaks = [[0] * 2 for i in range(len(prod1))]
    tminpeaks = [[0] * 2 for i in range(len(prod1))]
    inminpeaks = [[0] * 2 for i in range(len(prod1))]
    t1 = [[] for i in range(len(prod1))]
    t2 = [[] for i in range(len(prod1))]
    t3 = [[] for i in range(len(prod1))]
    t4 = [[] for i in range(len(prod1))]
    i1 = [[] for i in range(len(prod1))]
    i2 = [[] for i in range(len(prod1))]
    i3 = [[] for i in range(len(prod1))]
    i4 = [[] for i in range(len(prod1))]

    for i in range(len(prod1)):
        tmax[i] = max(t[i])
        tmin[i] = min(t[i])
        inmax[i] = np.multiply(t[i].index(tmax[i]), 2)
        inmin[i] = np.multiply(t[i].index(tmin[i]), 2)
        tempmax[i], _ = find_peaks(t[i])
        tempmin[i], _ = find_peaks(np.multiply(t[i], -1))
        inmaxpeaks[i] = np.multiply(tempmax[i], 2)
        inminpeaks[i] = np.multiply(tempmin[i], 2)
        for j in range(len(tempmax[i])):
            (tmaxpeaks[i])[j] = (t[i])[(tempmax[i])[j]]
            (tminpeaks[i])[j] = (t[i])[(tempmin[i])[j]]
    for i in range(len(prod1)):
        t1[i] = (tminpeaks[i])[0]
        t2[i] = (tmaxpeaks[i])[0]
        t3[i] = (tminpeaks[i])[1]
        t4[i] = (tmaxpeaks[i])[1]
        i1[i] = (inminpeaks[i])[0]
        i2[i] = (inmaxpeaks[i])[0]
        i3[i] = (inminpeaks[i])[1]
        i4[i] = (inmaxpeaks[i])[1]

    case = {"earheight": [(prod1[i])[0] for i in range(len(prod1))],
            "aslheight": [(prod1[i])[1] for i in range(len(prod1))],
            "torque": t,
            "inmax": inmax,
            "inmin": inmin,
            "tmax": tmax,
            "tmin": tmin,
            "inmaxpeaks": inmaxpeaks,
            "inminpeaks": inminpeaks,
            "tmaxpeaks": tmaxpeaks,
            "tminpeaks": tminpeaks,
            "t1": t1,
            "t2": t2,
            "t3": t3,
            "t4": t4,
            "i1": i1,
            "i2": i2,
            "i3": i3,
            "i4": i4}
    case = pd.DataFrame(case)
    case.to_pickle(ModelDir.DATA / "df_locked.pkl")

elif switch == 2:
    t = [[] for i in range(len(prod1))]
    a = 0
    b = 0
    while a < len(prod1):
        t[a] = [(torque["Torque"])[i] for i in range(b + 0, b + 91)]
        t[a] = [round((t[a])[i], 3) for i in range(len(range_c))]
        a = a + 1
        b = b + 91

    tmax = [[] for i in range(len(prod1))]
    tmin = [[] for i in range(len(prod1))]
    inmax = [[] for i in range(len(prod1))]
    inmin = [[] for i in range(len(prod1))]
    tempmax = [[] for i in range(len(prod1))]
    tempmin = [[] for i in range(len(prod1))]
    tmaxpeaks = [[0] * 2 for i in range(len(prod1))]
    inmaxpeaks = [[0] * 2 for i in range(len(prod1))]
    tminpeaks = [[0] * 2 for i in range(len(prod1))]
    inminpeaks = [[0] * 2 for i in range(len(prod1))]
    t1 = [[] for i in range(len(prod1))]
    t2 = [[] for i in range(len(prod1))]
    t3 = [[] for i in range(len(prod1))]
    t4 = [[] for i in range(len(prod1))]
    i1 = [[] for i in range(len(prod1))]
    i2 = [[] for i in range(len(prod1))]
    i3 = [[] for i in range(len(prod1))]
    i4 = [[] for i in range(len(prod1))]

    for i in range(len(prod1)):
        tmax[i] = max(t[i])
        tmin[i] = min(t[i])
        inmax[i] = np.multiply(t[i].index(tmax[i]), 2)
        inmin[i] = np.multiply(t[i].index(tmin[i]), 2)
        tempmax[i], _ = find_peaks(t[i])
        tempmin[i], _ = find_peaks(np.multiply(t[i], -1))
        inmaxpeaks[i] = np.multiply(tempmax[i], 2)
        inminpeaks[i] = np.multiply(tempmin[i], 2)
        for j in range(len(tempmax[i])):
            (tmaxpeaks[i])[j] = (t[i])[(tempmax[i])[j]]
            (tminpeaks[i])[j] = (t[i])[(tempmin[i])[j]]
    for i in range(len(prod1)):
        t1[i] = (tminpeaks[i])[0]
        t2[i] = (tmaxpeaks[i])[0]
        t3[i] = (tminpeaks[i])[1]
        t4[i] = (tmaxpeaks[i])[1]
        i1[i] = (inminpeaks[i])[0]
        i2[i] = (inmaxpeaks[i])[0]
        i3[i] = (inminpeaks[i])[1]
        i4[i] = (inmaxpeaks[i])[1]

    case = {"earheight": [(prod1[i])[0] for i in range(len(prod1))],
            "aslheight": [(prod1[i])[1] for i in range(len(prod1))],
            "torque": t,
            "inmax": inmax,
            "inmin": inmin,
            "tmax": tmax,
            "tmin": tmin,
            "inmaxpeaks": inmaxpeaks,
            "inminpeaks": inminpeaks,
            "tmaxpeaks": tmaxpeaks,
            "tminpeaks": tminpeaks,
            "t1": t1,
            "t2": t2,
            "t3": t3,
            "t4": t4,
            "i1": i1,
            "i2": i2,
            "i3": i3,
            "i4": i4}
    case = pd.DataFrame(case)
    case.to_pickle(ModelDir.DATA / "df_locked.pkl")

elif switch == 3:
    t = [[] for i in range(len(prod2))]
    a = 0
    b = 0
    while a < len(prod2):
        t[a] = (torque["Torque"])[b]
        if a in list(90*(i+1) for i in range(25)):
            b = b + 2822
            a = a + 1
        else:
            b = b + 1
            a = a + 1

    case = {"earheight": [(prod2[i])[0] for i in range(len(prod2))],
            "angle": [(prod2[i])[1] for i in range(len(prod2))],
            "torque": t}

    case = pd.DataFrame(case)
    #case.to_pickle(ModelDir.DATA / "df_locked.pkl")

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection='3d')
    a = 0
    b = 91
    for i in range(26):
        zdata = case["torque"].iloc[a+(i*91): b+(i*91)]
        xdata = case["angle"].iloc[a+(i*91): b+(i*91)]
        ydata = case["earheight"].iloc[a+(i*91): b+(i*91)]
        ax.scatter3D(xdata, ydata, zdata)
    ax.set_xlabel('Angle [deg]', fontsize=10)
    ax.set_ylabel('Parameter A [mm]', fontsize=10)
    ax.set_zlabel('Torque [Nm]', fontsize=10)
    ax.minorticks_on()
    ax.set_ylim(3, 0)
    ax.set_ylim(3, 0)
    ax.tick_params(labelsize=10)
    ax.view_init(elev=20, azim=250)
    #plt.savefig(ModelDir.MEDIA / "i4_locked3d.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 4:
    t = [[] for i in range(len(prod3))]
    a = 0
    b = 45136
    while a < len(prod3):
        t[a] = (torque["Torque"])[b]
        b = b + 1
        a = a + 1
    print(b)
    case = {"aslheight": [(prod3[i])[0] for i in range(len(prod3))],
            "angle": [(prod3[i])[1] for i in range(len(prod3))],
            "torque": t}

    case = pd.DataFrame(case)
    #case.to_pickle(ModelDir.DATA / "df_locked.pkl")

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection='3d')
    a = 55
    b = 65
    for i in range(len(range_b)):
        zdata = case["torque"].iloc[a+(i*91): b+(i*91)]
        xdata = case["angle"].iloc[a+(i*91): b+(i*91)]
        ydata = case["aslheight"].iloc[a+(i*91): b+(i*91)]
        ax.scatter3D(xdata, ydata, zdata)
    ax.set_xlabel('Angle [deg]', fontsize=10)
    ax.set_ylabel('Parameter C [mm]', fontsize=10)
    ax.set_zlabel('Torque [Nm]', fontsize=10)
    ax.minorticks_on()
    ax.set_ylim(3, 0)
    ax.set_ylim(3, 0)
    ax.tick_params(labelsize=10)
    ax.view_init(elev=20, azim=250)
    #plt.savefig(ModelDir.MEDIA / "i4_locked3d.png", bbox_inches="tight", dpi=650)
    plt.show()

else:
    pass