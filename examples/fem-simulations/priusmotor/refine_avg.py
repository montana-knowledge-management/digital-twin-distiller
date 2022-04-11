import numpy as np
from numpy import linspace
from itertools import product
import json
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from digital_twin_distiller import ModelDir

ModelDir.set_base(__file__)

switch = 0
if switch == 0:

    range_b0 = 0
    range_b1 = -60
    nsteps_b = 61

    range_c0 = 100
    range_c1 = 250
    nsteps_c = 4

    iterlist = [[], [], []]
    for i in range(nsteps_c):
        for j in range(451*4):
            range_a0 = 0 + j / 10
            range_a1 = 15 + j / 10
            nsteps_a = nsteps_b

            range_a = linspace(range_a0, range_a1, nsteps_a)
            range_b = linspace(range_b0, range_b1, nsteps_b)
            range_c = linspace(range_c0, range_c1, nsteps_c)

            for k in range(nsteps_b):
                iterlist[0].append(round(range_a[k], 3))
                iterlist[1].append(round(range_b[k], 3))
                iterlist[2].append(round(range_c[i], 3))

    f = open(ModelDir.DATA / f'avg.json')
    avg = json.load(f)

    t = [[] for i in range(451*4)]
    a = 0
    b = 0
    while a < 451*4:
        t[a] = [(avg["Torque"])[i] for i in range(b + 0, b + 61)]
        a = a + 1
        b = b + 61

    tavg = [[] for i in range(451*4)]
    tmax = [[] for i in range(451 * 4)]
    tmin = [[] for i in range(451 * 4)]
    twav = [[] for i in range(451 * 4)]
    rota = [[] for i in range(451)]
    temp = np.multiply((range(451*4)), 0.4)
    for i in range(451*4):
        tavg[i] = np.mean(t[i])
        tmax[i] = max(t[i])
        tmin[i] = min(t[i])
        twav[i] = (tmax[i] - tavg[i]) / tavg[i]
    for i in range(451):
        rota[i] = temp[i]
    for i in range(3):
        for j in range(451):
            rota.append(rota[j])
    for i in range(451*4):
        rota[i] = round(rota[i], 3)
        tavg[i] = round(tavg[i], 3)
        tmax[i] = round(tmax[i], 3)
        tmin[i] = round(tmin[i], 3)
        twav[i] = round(twav[i], 3)

    res = {"rotorangle": rota,
           "tavg": tavg,
           "tmax": tmax,
           "tmin": tmin,
           "twav": twav}

    res = pd.DataFrame(res)
    res.to_pickle(ModelDir.DATA / "df_avg.pkl")
    print(res)