import numpy as np
from numpy import linspace
from itertools import product
import json
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from digital_twin_distiller import ModelDir

ModelDir.set_base(__file__)

switch = 3
if switch == 0:

    range_a0 = 0
    range_a1 = 250
    nsteps_a = 251

    range_b0 = 0
    range_b1 = 45
    nsteps_b = 91

    range_a = linspace(range_a0, range_a1, nsteps_a)
    range_b = linspace(range_b0, range_b1, nsteps_b)
    prod = list(product(range_a, range_b))

    f = open(ModelDir.DATA / f'locked_i.json')
    rotate = json.load(f)

    t = [[] for i in range(len(range_a))]
    a = 0
    b = 0
    while a < len(range_a):
        t[a] = [(rotate["Torque"])[i] for i in range(b + 0, b + 91)]
        a = a + 1
        b = b + 91

    tempmax = [[] for i in range(len(range_a))]
    tmaxpeaks = [[0] * 10 for i in range(len(range_a))]
    inmaxpeaks = [[0] * 10 for i in range(len(range_a))]
    tempmin = [[] for i in range(len(range_a))]
    tminpeaks = [[0] * 10 for i in range(len(range_a))]
    inminpeaks = [[0] * 10 for i in range(len(range_a))]
    t1 = [[] for i in range(len(range_a))]
    t2 = [[] for i in range(len(range_a))]
    t3 = [[] for i in range(len(range_a))]
    t01 = [[] for i in range(len(range_a))]
    t02 = [[] for i in range(len(range_a))]
    t03 = [[] for i in range(len(range_a))]
    tav = [[] for i in range(len(range_a))]

    for i in range(len(range_a)):
        tempmax[i], _ = find_peaks(t[i])
        tempmin[i], _ = find_peaks(np.multiply(t[i], -1))
        inmaxpeaks[i] = np.multiply(tempmax[i], 2)
        inminpeaks[i] = np.multiply(tempmin[i], 2)
        for j in range(len(tempmax[i])):
            (tmaxpeaks[i])[j] = (t[i])[(tempmax[i])[j]]
        for j in range(len(tempmin[i])):
            (tminpeaks[i])[j] = (t[i])[(tempmin[i])[j]]
        t1[i] = (tmaxpeaks[i])[0]
        t2[i] = (tminpeaks[i])[1]
        t3[i] = (tmaxpeaks[i])[1]
        tav[i] = np.multiply(t1[i] + t2[i] + t3[i], (1 / 3))
        t01[i] = t1[i] - tav[i]
        t02[i] = t2[i] - tav[i]
        t03[i] = t3[i] - tav[i]

    print(t[250])

    res = {"current": range_a,
           "rotorangle": [np.multiply(range_b, 4) * (i+1)/(i+1) for i in range(len(range_a))],
           "torque": [t[a] for a in range(len(range_a))],
           "t1peaks": tmaxpeaks,
           "i1peaks": inmaxpeaks,
           "t2peaks": tminpeaks,
           "i2peaks": inminpeaks,
           "t1": t1,
           "t2": t2,
           "t3": t3,
           "t01": t01,
           "t02": t02,
           "t03": t03,
           "tav": tav}

    res = pd.DataFrame(res)
    res.to_pickle(ModelDir.DATA / "df_rotate0.pkl")

elif switch == 1:

    range_a0 = 50
    range_a1 = 250
    nsteps_a = 6

    range_b0 = 30.00
    range_b1 = 45.00
    nsteps_b = 61

    range_c0 = 0
    range_c1 = -60
    nsteps_c = 61

    range_a = linspace(range_a0, range_a1, nsteps_a)
    range_b = linspace(range_b0, range_b1, nsteps_b)
    range_c = linspace(range_c0, range_c1, nsteps_c)
    prod = list(product(range_a, range_b))

    f = open(ModelDir.DATA / f'rotate_t2.json')
    rotate = json.load(f)

    t = [[] for i in range(len(range_a))]
    tav = [[] for i in range(len(range_a))]
    tmax = [[] for i in range(len(range_a))]
    twav = [[] for i in range(len(range_a))]
    a = 0
    b = 0
    while a < len(range_a):
        t[a] = [(rotate["Torque"])[i] for i in range(b + 0, b + 61)]
        tav[a] = np.average(t[a])
        tmax[a] = max(t[a])
        twav[a] = tmax[a] - tav[a]
        a = a + 1
        b = b + 61

    res = {"current": range_a,
           "rotorangle": [np.multiply(range_c, -1) * (i + 1) / (i + 1) for i in range(len(range_a))],
           "torque": [t[a] for a in range(len(range_a))],
           "tav": [tav[a] for a in range(len(range_a))],
           "tmax": tmax,
           "twav": twav}
    res = pd.DataFrame(res)
    res.to_pickle(ModelDir.DATA / "df_rotate_t2.pkl")

    plt.scatter(res["current"], res["tmax"])
    plt.show()

elif switch == 2:

    range_a0 = 0
    range_a1 = 250
    nsteps_a = 251

    range_b0 = 30.00
    range_b1 = 45.00
    nsteps_b = 61

    range_c0 = 0
    range_c1 = -60
    nsteps_c = 61

    range_a = linspace(range_a0, range_a1, nsteps_a)
    range_b = linspace(range_b0, range_b1, nsteps_b)
    range_c = linspace(range_c0, range_c1, nsteps_c)
    prod = list(product(range_a, range_b))

    f = open(ModelDir.DATA / f'rotateit2.json')
    rotate = json.load(f)

    t = [[] for i in range(len(range_a))]
    tav = [[] for i in range(len(range_a))]
    tmax = [[] for i in range(len(range_a))]
    twav = [[] for i in range(len(range_a))]
    a = 0
    b = 0
    while a < len(range_a):
        t[a] = [(rotate["Torque"])[i] for i in range(b + 0, b + 61)]
        tav[a] = np.mean(t[a])
        tmax[a] = max(t[a])
        twav[a] = tmax[a] - tav[a]
        a = a + 1
        b = b + 61

    res = {"current": range_a,
           "rotorangle": [np.multiply(range_c, -1) * (i + 1) / (i + 1) for i in range(len(range_a))],
           "torque": [t[a] for a in range(len(range_a))],
           "tav": [tav[a] for a in range(len(range_a))],
           "tmax": tmax,
           "twav": twav}
    res = pd.DataFrame(res)
    res.to_pickle(ModelDir.DATA / "df_rotateit2.pkl")

    f = open(ModelDir.DATA / f'rotateit3.json')
    rotate = json.load(f)

    t = [[] for i in range(len(range_a))]
    tav = [[] for i in range(len(range_a))]
    tmax = [[] for i in range(len(range_a))]
    twav = [[] for i in range(len(range_a))]
    a = 0
    b = 0
    while a < len(range_a):
        t[a] = [(rotate["Torque"])[i] for i in range(b + 0, b + 61)]
        tav[a] = np.mean(t[a])
        tmax[a] = max(t[a])
        twav[a] = tmax[a] - tav[a]
        a = a + 1
        b = b + 61

    res = {"current": range_a,
           "rotorangle": [np.multiply(range_c, -1) * (i + 1) / (i + 1) for i in range(len(range_a))],
           "torque": [t[a] for a in range(len(range_a))],
           "tav": [tav[a] for a in range(len(range_a))],
           "tmax": tmax,
           "twav": twav}
    res = pd.DataFrame(res)
    res.to_pickle(ModelDir.DATA / "df_rotateit3.pkl")

elif switch == 3:

    f = open(ModelDir.DATA / f'testrot.json')
    rotate = json.load(f)

    t = [[] for i in range(11)]
    tav = [[] for i in range(11)]
    a = 0
    b = 0
    while a < 11:
        t[a] = [(rotate["Torque"])[i] for i in range(b + 0, b + 16)]
        tav[a] = np.mean(t[a])
        a = a + 1
        b = b + 16
    range_c = linspace(0, 60, 16)
    res = {"rotorangle": [np.multiply(range_c, 1) * (i + 1) / (i + 1) for i in range(11)],
           "torque": [t[a] for a in range(11)]}
    res = pd.DataFrame(res)
    a0 = 0
    a1 = 5
    b = 16
    fig = plt.subplots(figsize=(6, 4))
    for c in range(a0, a1):
        plt.plot([((res["rotorangle"])[c])[d] for d in range(b)], [((res["torque"])[c])[d] for d in range(b)],
                 label=str(c))
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "x.png", bbox_inches="tight", dpi=650)
    plt.show()