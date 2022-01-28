from itertools import product as pr
from numpy import linspace

range_a0 = 0
range_a1 = 1
nsteps_a = 2

range_b0 = 0
range_b1 = 1
nsteps_b = 2

range_c0 = 0
range_c1 = 1
nsteps_c = 2

range_a = linspace(range_a0, range_a1, nsteps_a)
range_b = linspace(range_b0, range_b1, nsteps_a)
range_c = linspace(range_c0, range_c1, nsteps_c)

prod = list(pr(range_a, range_b, range_c))

for i in range(len(prod)):
    earheight = (prod[i])[0]
    aslheight = (prod[i])[1]
    alpha = (prod[i])[2]

