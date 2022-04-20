from statistics import mean

import numpy as np
import scipy.signal
from numpy import linspace
from itertools import product
import json
import pandas as pd
from digital_twin_distiller import ModelDir
import matplotlib.pyplot as plt

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

f = open(ModelDir.DATA / f'cogging_torque.json')
torque = json.load(f)

res = {}
res["earheight"] = [(prod[i])[0] for i in range(len(prod))]
res["aslheight"] = [(prod[i])[1] for i in range(len(prod))]
res["rotorangle"] = [(prod[i])[2] for i in range(len(prod))]
res["torque"] = [(torque["Torque"])[i] for i in range(len(prod))]
res = pd.DataFrame(res)

switch = 4
if switch == 0:
    t = [[] for i in range(len(prod1))]
    a = 0
    b = 0
    while a < len(prod1):
        t[a] = [(torque["Torque"])[i] for i in range(b+0, b+76)]
        a = a + 1
        b = b + 76

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
                (inpeak[a])[b] = round(b / 75 * 7.5, 1)
        inpeak[a] = [i for i in inpeak[a] if i != 0]
        tpeak[a] = [i for i in tpeak[a] if i != 0]

    case = {"earheight": [(prod1[i])[0] for i in range(len(prod1))],
            "aslheight": [(prod1[i])[1] for i in range(len(prod1))],
            "coggingtorque": t,
            "torquepeak": tpeak,
            "peakindex": inpeak}
    case = pd.DataFrame(case)
    print(case["torquepeak"])

    case.to_pickle(ModelDir.DATA / "df_cogging.pkl")

elif switch == 1:
    t = [[] for i in range(nsteps_a)]
    a = 0
    b = 0
    while a < nsteps_a:
        t[a] = [(torque["Torque"])[i] for i in range(b + 0, b + 76)]
        a = a + 1
        b = b + 2356
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
                (inpeak[a])[b] = round(b / 75 * 7.5, 1)
        inpeak[a] = [i for i in inpeak[a] if i != 0]
        tpeak[a] = [i for i in tpeak[a] if i != 0]

    case = {"earheight": range_a,
            "coggingtorque": t,
            "torquepeak": tpeak,
            "peakindex": inpeak}
    case = pd.DataFrame(case)
    print(case)

    case.to_pickle(ModelDir.DATA / "df_cogging.pkl")

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
    tdelta3 = [[] for i in range(len(prod1))]
    for i in range(len(prod1)):
        tmax[i] = max(t[i])
        tmin[i] = min(t[i])
        inmax[i] = np.multiply(t[i].index(tmax[i]), 0.2) +15
        inmin[i] = np.multiply(t[i].index(tmin[i]), 0.2) +15
        tdelta1[i] = abs(tmax[i] - tmin[i])
        tdelta2[i] = tmax[i] - abs(tmin[i])
        tdelta3[i] = mean(t[i])


    case = {"earheight": [(prod1[i])[0] for i in range(len(prod1))],
            "aslheight": [(prod1[i])[1] for i in range(len(prod1))],
            "coggingtorque": t,
            "inmaxpeak": inmax,
            "inminpeak": inmin,
            "tmaxpeak": tmax,
            "tminpeak": tmin,
            "tdelta1": tdelta1,
            "tdelta2": tdelta2,
            "tdelta3": tdelta3
            }
    case = pd.DataFrame(case)
    case['inminpeak'] = case["inminpeak"].replace({15: np.nan})
    case['inminpeak'] = case["inminpeak"].replace({30: np.nan})
    case["tminpeak"].values[case['tminpeak'] > -0.05] = np.nan

    print(case["inminpeak"])

    case.to_pickle(ModelDir.DATA / "df_cogging.pkl")
else:
    pass