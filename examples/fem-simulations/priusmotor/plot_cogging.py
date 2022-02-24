from itertools import product

import matplotlib.pyplot as plt
from digital_twin_distiller import ModelDir
from numpy import linspace
import pandas as pd
from ast import literal_eval

ModelDir.set_base(__file__)

range_a0 = 0.5
range_a1 = 3.0
nsteps_a = 26

range_b0 = 0.0
range_b1 = 3.0
nsteps_b = 31

range_c0 = 0.0
range_c1 = 7.5
nsteps_c = 76

range_a = linspace(range_a0, range_a1, nsteps_a)
range_b = linspace(range_b0, range_b1, nsteps_b)
range_c = linspace(range_c0, range_c1, nsteps_c)

prod = list(product(range_a, range_b, range_c))
range_prod = linspace (0, len(prod), len(prod)+1)
prod1 = list(product(range_a, range_b))

range_c = linspace(range_c0, range_c1, nsteps_c)

case = pd.read_pickle(ModelDir.DATA / "df_cogging.pkl")

#fig, ax1 = plt.subplots()
#plt.title("Peak of the cogging torque and the rotor angle")
#ax2 = ax1.twinx()
#for xe, ye in zip(range_a, case["torquepeak"]):
    #ax1.scatter([xe]*len(ye), ye,)
#for xe, ye in zip(range_a, case["peakindex"]):
    #ax2.scatter([xe]*len(ye), ye)
#ax1.set_xlabel('Variable parameter [mm]', color='r')
#ax1.set_ylabel('Peak cogging torque [Nm]', color='g')
#ax2.set_ylabel('Rotor angle [°]', color='b')
#plt.show()


for xe, ye in zip(range(len(prod1)), case["torquepeak"]):
    plt.scatter([xe]*len(ye), ye)
plt.show()
for xe, ye in zip(range(len(prod1)), case["peakindex"]):
    plt.scatter([xe]*len(ye), ye)
plt.show()

a = 0
b = 26
for i in range(a, b):
    plt.plot(range_c, (case["coggingtorque"])[i], lw=2)
plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
plt.minorticks_on()
plt.xlabel("rotorangle [deg]")
plt.ylabel("Torque [Nm]")
plt.show()


