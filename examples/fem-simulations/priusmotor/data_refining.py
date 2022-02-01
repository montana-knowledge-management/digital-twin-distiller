from numpy import linspace
from itertools import product
import json

from digital_twin_distiller import ModelDir

ModelDir.set_base(__file__)

range_a0 = 0.5
range_a1 = 3.5
nsteps_a = 31

range_b0 = 0.0
range_b1 = 3.0
nsteps_b = 31

range_c0 = 0.0
range_c1 = 7.5
nsteps_c = 16

range_a = linspace(range_a0, range_a1, nsteps_a)
range_b = linspace(range_b0, range_b1, nsteps_b)
range_c = linspace(range_c0, range_c1, nsteps_c)

prod = list(product(range_a, range_b, range_c))

f = open(ModelDir.DATA / f'cogging_torque.json')
torque = json.load(f)

