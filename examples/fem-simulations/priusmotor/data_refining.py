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
range_c1 = 7.5
nsteps_c = 16

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

t = [[] for i in range(961)]
a = 0
b = 0
while a <= 960:
    t[a] = [(torque["Torque"])[i] for i in range(b+0, b+16)]
    a = a + 1
    b = b + 16

case = {}
case["earheight"] = [(prod1[i])[0] for i in range(len(prod1))]
case["aslheight"] = [(prod1[i])[1] for i in range(len(prod1))]
case["torque"] = t
case = pd.DataFrame(case)

null = []
comp = [-0.0 for i in range(931)]
null = [((case["torque"])[i])[0] for i in range(931)]
y = []
for l1,l2 in zip(comp,null):
    if l1 == l2:
        x = 1
    else:
        x = 0
    y.append(x)
y = pd.DataFrame(y)

z = y.loc[y[0] == 1]
print(z)

print(case.loc[24])

#for i in range(931):
    #plt.plot(range_c, (case["torque"])[i], lw=2)
#plt.plot(range_c, (case["torque"])[1], lw=2)
#plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
#plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
#plt.minorticks_on()
#plt.xlabel("rotorangle [deg]")
#plt.ylabel("Cogging Torque [Nm]")
#plt.show()