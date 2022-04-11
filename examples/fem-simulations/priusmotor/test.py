import numpy as np
from numpy import linspace

range_a0 = 100
range_a1 = 250
nsteps_a = 4

range_b0 = 0.0
range_b1 = 3.0
nsteps_b = 31

range_c0 = 0.0
range_c1 = -60.0
nsteps_c = 61

iterlist = [[], [], [], []]
range_a = linspace(range_a0, range_a1, nsteps_a)
range_b = linspace(range_b0, range_b1, nsteps_b)
range_c = linspace(range_c0, range_c1, nsteps_c)

for a in range(nsteps_a):
    for b in range(nsteps_b):
        for i in range(46):
            range_d0 = 0 + i
            range_d1 = 15 + i
            nsteps_d = nsteps_c

            range_d = linspace(range_d0, range_d1, nsteps_d)

            for j in range(nsteps_d):
                iterlist[0].append(round(range_d[j], 3))
                iterlist[1].append(round(range_c[j], 3))
                iterlist[2].append(round(range_b[b], 3))
                iterlist[3].append(round(range_a[a], 3))



