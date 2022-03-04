from itertools import product

import matplotlib.pyplot as plt
from digital_twin_distiller import ModelDir
from numpy import linspace
import pandas as pd

ModelDir.set_base(__file__)

range_a0 = 0.5
range_a1 = 3.0
nsteps_a = 26

range_b0 = 0.0
range_b1 = 3.0
nsteps_b = 31

range_c0 = 0.0
range_c1 = 45
nsteps_c = 91

range_a = linspace(range_a0, range_a1, nsteps_a)
range_b = linspace(range_b0, range_b1, nsteps_b)
range_c = linspace(range_c0, range_c1, nsteps_c)

prod = list(product(range_a, range_b, range_c))
range_prod = linspace (0, len(prod), len(prod)+1)
prod1 = list(product(range_a, range_b))

range_c = linspace(range_c0, range_c1, nsteps_c)

case = pd.read_pickle(ModelDir.DATA / "df_locked.pkl")

switch = 1
if switch == 1:
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(projection='3d')
    c = 0
    a = 0
    b = 806
    for i in range(a, b):
        zdata = (case["tmaxpeak"].loc[range(a, b)]).tolist()
        xdata = (case["earheight"].loc[range(a, b)]).tolist()
        ydata = (case["aslheight"].loc[range(a, b)]).tolist()
        ax.scatter3D(xdata, ydata, zdata)
        ax.set_xlabel('Parameter A [mm]', fontsize=10)
        ax.set_ylabel('Parameter C [mm]', fontsize=10)
        ax.set_zlabel('Torque [Nm]', fontsize=10)
        ax.minorticks_on()
        ax.tick_params(labelsize=10)
        plt.gca().invert_yaxis()
        ax.view_init(elev=20., azim=140)
    #plt.savefig(ModelDir.MEDIA / "cogging3d.png", bbox_inches="tight", dpi=650)
    plt.show()