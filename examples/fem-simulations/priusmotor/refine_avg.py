import numpy as np
import pandas as pd
from numpy import linspace
import json
import matplotlib.pyplot as plt


from digital_twin_distiller import ModelDir

ModelDir.set_base(__file__)

switch = 0

if switch == 0:

    range_a0 = 0.5
    range_a1 = 3.0
    nsteps_a = 26

    range_b0 = 0.0
    range_b1 = -60.0
    nsteps_b = 61

    iterlist = [[], [], []]
    range_a = linspace(range_a0, range_a1, nsteps_a)
    range_b = linspace(range_b0, range_b1, nsteps_b)

    for a in range(nsteps_a):
        for i in range(46):
            range_c0 = 0 + i
            range_c1 = 15 + i
            nsteps_c = nsteps_b

            range_c = linspace(range_c0, range_c1, nsteps_c)

            for j in range(nsteps_c):
                iterlist[0].append(round(range_c[j], 3))
                iterlist[1].append(round(range_b[j], 3))
                iterlist[2].append(round(range_a[a], 3))

    f = open(ModelDir.DATA / f'avgear_100.json')
    avg = json.load(f)

    ran = 46 * nsteps_a

    t = [[] for i in range(ran)]
    a = 0
    b = 0
    while a < (ran):
        t[a] = [(avg["Torque"])[i] for i in range(b + 0, b + 61)]
        a = a + 1
        b = b + 61

    tavg = [[] for i in range(ran)]
    #tmax = [[] for i in range(ran)]
    #tmin = [[] for i in range(ran)]
    #twav = [[] for i in range(ran)]
    rota = [[] for i in range(46)]
    for i in range(ran):
        tavg[i] = np.mean(t[i])
    #print(len(tavg))

        #tmax[i] = max(t[i])
        #tmin[i] = min(t[i])
        #twav[i] = (tmax[i] - tavg[i]) / tavg[i]
    for i in range(46):
        rota[i] = 4*i
    for i in range(25):
        for j in range(46):
            rota.append(rota[j])
    #print(len(rota))
    for i in range(ran):
        rota[i] = round(rota[i], 3)
        tavg[i] = round(tavg[i], 3)
        #tmax[i] = round(tmax[i], 3)
        #tmin[i] = round(tmin[i], 3)
        #twav[i] = round(twav[i], 3)

    x = [[] for i in range(nsteps_a)]
    y = [[] for i in range(nsteps_a)]
    for i in range(26):
        for j in range(46):
            x[i].append(rota[j])
    a=0
    b=0
    while a < 26:
        for c in range(46):
            y[a].append(tavg[c+b])
        a=a+1
        b=b+46

    res = {"rotorangle": x,
          "tavg": y}

    res = pd.DataFrame(res)
    res.to_pickle(ModelDir.DATA / "df_avgear100.pkl")
    print(res)