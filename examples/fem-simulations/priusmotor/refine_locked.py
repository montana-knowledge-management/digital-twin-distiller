from numpy import linspace
from itertools import product
import json
import pandas as pd
import matplotlib.pyplot as plt

from digital_twin_distiller import ModelDir, setup_matplotlib

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

f = open(ModelDir.DATA / f'locked_rotor.json')
locked = json.load(f)

res = {"earheight": [(prod[i])[0] for i in range(len(prod))],
       "aslheight": [(prod[i])[1] for i in range(len(prod))],
       "rotorangle": [(prod[i])[2] for i in range(len(prod))],
       "torque": [(locked["Torque"])[i] for i in range(len(prod))]}
res = pd.DataFrame(res)

switch = 1
if switch == 0:
    t = [[] for i in range(len(prod1))]
    a = 0
    b = 0
    while a < len(prod1):
        t[a] = [(locked["Torque"])[i] for i in range(b+0, b+25)]
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
                (inpeak[a])[b] = b * 1.875
        inpeak[a] = [i for i in inpeak[a] if i != 0]
        tpeak[a] = [i for i in tpeak[a] if i != 0]

    case = {"earheight": [(prod1[i])[0] for i in range(len(prod1))],
            "aslheight": [(prod1[i])[1] for i in range(len(prod1))],
            "torque": t,
            "torquepeak": tpeak,
            "peakindex": inpeak}
    case = pd.DataFrame(case)
    print(case)

    case.to_pickle(ModelDir.DATA / "df_locked.pkl")
else:
    t = [[] for i in range(nsteps_a)]
    a = 0
    b = 0
    while a < nsteps_a:
        t[a] = [(locked["Torque"])[i] for i in range(b + 0, b + 25)]
        a = a + 1
        b = b + 775
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
                (inpeak[a])[b] = b * 1.875
        inpeak[a] = [i for i in inpeak[a] if i != 0]
        tpeak[a] = [i for i in tpeak[a] if i != 0]

    case = {"earheight": range_a,
            "torque": t,
            "torquepeak": tpeak,
            "peakindex": inpeak}
    case = pd.DataFrame(case)
    print(case)

    case.to_pickle(ModelDir.DATA / "df_locked.pkl")