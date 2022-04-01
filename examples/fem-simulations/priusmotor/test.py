import math
from itertools import product as pr
from numpy import linspace

from digital_twin_distiller import Node

range_b0 = 0
range_b1 = -60
nsteps_b = 61

range_c0 = 100
range_c1 = 250
nsteps_c = 4

iterlist = [[], [], []]
for i in range(nsteps_c):
    for j in range(451):
        range_a0 = 0 + j/10
        range_a1 = 15 + j/10
        nsteps_a = nsteps_b

        range_a = linspace(range_a0, range_a1, nsteps_a)
        range_b = linspace(range_b0, range_b1, nsteps_b)
        range_c = linspace(range_c0, range_c1, nsteps_c)

        for k in range(nsteps_b):
            iterlist[0].append(round(range_a[k], 3))
            iterlist[1].append(round(range_b[k], 3))
            iterlist[2].append(round(range_c[i], 3))


print(iterlist)

