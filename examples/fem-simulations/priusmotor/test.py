import numpy as np
from numpy import linspace

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

print(len(iterlist[0]))

