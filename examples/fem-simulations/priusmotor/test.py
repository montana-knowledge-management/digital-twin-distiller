import math
from itertools import product as pr
from numpy import linspace

from digital_twin_distiller import Node

range_a0 = 1
range_a1 = 2
nsteps_a = 2

range_b0 = 1
range_b1 = 2
nsteps_b = 2

range_c0 = 1
range_c1 = 2
nsteps_c = 2

range_a = linspace(range_a0, range_a1, nsteps_a)
range_b = linspace(range_b0, range_b1, nsteps_a)
range_c = linspace(range_c0, range_c1, nsteps_c)

prod = list(pr(range_a, range_b, range_c))
iterat = list(range(len(prod)))

def pol2cart(rho: float, phi: float):
    x = rho * math.cos(math.radians(phi))
    y = rho * math.sin(math.radians(phi))
    return x, y

temp1 = pol2cart(1, 1)
print(temp1)