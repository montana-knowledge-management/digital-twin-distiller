from itertools import product
import matplotlib.pyplot as plt
import numpy as np
from digital_twin_distiller import ModelDir
from numpy import linspace
import pandas as pd
from matplotlib.lines import Line2D
import numpy as np

ModelDir.set_base(__file__)

switch = 1
if switch == 0:
    res = pd.read_pickle(ModelDir.DATA / "df_rotate0.pkl")
    a = 1
    b = 91
    fig = plt.subplots(figsize=(6, 4))
    for c in range(a):
        plt.plot([((res["rotorangle"])[c])[d] for d in range(b)], [((res["torque"])[c])[d] for d in range(b)], label="Cogging torque")
    plt.xlabel('Electrical angle [deg]', fontsize=10)
    plt.ylabel('Torque [Nm]', fontsize=10)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(0, 200, step=20))
    plt.yticks(fontsize=10)
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "PEMC_T1.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 1:
    res = pd.read_pickle(ModelDir.DATA / "df_rotate0.pkl")

    a0 = 100
    a1 = 265
    a2 = 15
    b = 91
    fig = plt.figure(figsize=(6, 4))
    for c,g in zip(range(a0, a1, a2), range(100, 265, 15)):
        e = len(((res["i1peaks"])[c]))
        f = len(((res["i2peaks"])[c]))
        plt.plot([((res["rotorangle"])[c])[d] for d in range(b)], [((res["torque"])[c])[d] for d in range(b)], label=str(g) +"A")
        plt.scatter([((res["i1peaks"])[c])[d] for d in range(e)], [((res["t1peaks"])[c])[d] for d in range(e)], c='r')
        plt.scatter(((res["i2peaks"])[c])[1], ((res["t2peaks"])[c])[1], c='b')
    plt.xlabel('Electrical angle [deg]', fontsize=10)
    plt.ylabel('Torque [Nm]', fontsize=10)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(0, 200, step=20))
    plt.yticks(fontsize=10)
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "PEMC_T2.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 2:
    res = pd.read_pickle(ModelDir.DATA / "df_rotate0.pkl")

    plt.scatter((res["current"])[100:250], (res["t01"])[100:250], label = "T31")
    plt.scatter((res["current"])[100:250], (res["t02"])[100:250], label = "T32")
    plt.scatter((res["current"])[100:250], (res["t03"])[100:250], label = "T12")
    plt.xlabel('Current [A]', fontsize=10)
    plt.ylabel('Torque [Nm]', fontsize=10)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(100, 265, step=15))
    plt.yticks(fontsize=10)
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "PEMC_T3.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 3:
    fig = plt.figure(figsize=(6, 4))
    res1 = pd.read_pickle(ModelDir.DATA / "df_rotateit2.pkl")
    res2 = pd.read_pickle(ModelDir.DATA / "df_rotateit3.pkl")
    res3 = pd.read_pickle(ModelDir.DATA / "df_rotateit2.pkl")
    plt.plot(res1["current"], res1["tav"], label= r'$ \beta $' + "=" + str(120) + "deg")
    plt.plot(res2["current"], res2["tav"], label= r'$ \beta $' + "=" + str(132) + "deg")
    plt.plot(res3["current"], res3["tav"], label= r'$ \beta $' + "=" + str(148) + "deg")
    plt.xlabel('Current [A]', fontsize=10)
    plt.ylabel('Torque [Nm]', fontsize=10)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "PEMC_T4.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 4:
    fig = plt.figure(figsize=(6, 4))
    res1 = pd.read_pickle(ModelDir.DATA / "df_rotateit2.pkl")
    res2 = pd.read_pickle(ModelDir.DATA / "df_rotateit3.pkl")
    res3 = pd.read_pickle(ModelDir.DATA / "df_rotateit2.pkl")
    plt.plot(res1["current"], res1["twav"], label= r'$ \beta $' + "=" + str(120) + "deg")
    plt.plot(res2["current"], res2["twav"], label= r'$ \beta $' + "=" + str(132) + "deg")
    plt.plot(res3["current"], res3["twav"], label= r'$ \beta $' + "=" + str(148) + "deg")
    plt.xlabel('Current [A]', fontsize=10)
    plt.ylabel('Torque [Nm]', fontsize=10)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "PEMC_T5.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 5:
    res = pd.read_pickle(ModelDir.DATA / "df_rotate0.pkl")
    with open(ModelDir.DATA / 'locked250.csv', 'r', encoding='utf-8') as f:
        res_ = pd.read_csv(f)
        alpha250 = res_["x"]
        T250 = res_["y"]
    with open(ModelDir.DATA / 'locked200.csv', 'r', encoding='utf-8') as f:
        res_ = pd.read_csv(f)
        alpha200 = res_["x"]
        T200 = res_["y"]
    with open(ModelDir.DATA / 'locked150.csv', 'r', encoding='utf-8') as f:
        res_ = pd.read_csv(f)
        alpha150 = res_["x"]
        T150 = res_["y"]
    with open(ModelDir.DATA / 'locked100.csv', 'r', encoding='utf-8') as f:
        res_ = pd.read_csv(f)
        alpha100 = res_["x"]
        T100 = res_["y"]

    a0 = 100
    a1 = 300
    a2 = 50
    b = 91
    fig = plt.figure(figsize=(6, 4))
    for c, g in zip(range(a0, a1, a2), range(100, 300, 50)):
        plt.plot([((res["rotorangle"])[c])[d] for d in range(b)], [((res["torque"])[c])[d] for d in range(b)],
                 label=str(g) + "A")
        plt.scatter(alpha100, T100, c = "b")
        plt.scatter(alpha150, T150, c = "orange")
        plt.scatter(alpha200, T200, c = "g")
        plt.scatter(alpha250, T250, c = "r")
    plt.xlabel('Electrical angle [deg]', fontsize=10)
    plt.ylabel('Torque [Nm]', fontsize=10)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(0, 200, step=20))
    plt.yticks(fontsize=10)
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "PEMC_T6.png", bbox_inches="tight", dpi=650)
    plt.show()

if switch == 6:
    res1 = pd.read_pickle(ModelDir.DATA / "df_rotateit2.pkl")
    res2 = pd.read_pickle(ModelDir.DATA / "df_rotateit3.pkl")
    a0 = 200
    a1 = 251
    b = 61
    fig = plt.subplots(figsize=(6, 4))
    for c in range(a0, a1):
        plt.plot([((res1["rotorangle"])[c])[d] for d in range(b)], [((res1["torque"])[c])[d] for d in range(b)],
                 label="1", c="r")
        plt.plot([((res2["rotorangle"])[c])[d] for d in range(b)], [((res2["torque"])[c])[d] for d in range(b)],
                 label="2", c ="b")
    plt.xlabel('Electrical angle [deg]', fontsize=10)
    plt.ylabel('Torque [Nm]', fontsize=10)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(0, 80, step=20))
    plt.yticks(fontsize=10)
    #plt.legend()
    plt.savefig(ModelDir.MEDIA / "PEMC_T1.png", bbox_inches="tight", dpi=650)
    plt.show()
