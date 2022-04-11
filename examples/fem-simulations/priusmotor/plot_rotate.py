from itertools import product
import matplotlib.pyplot as plt
import numpy as np
from digital_twin_distiller import ModelDir
from numpy import linspace
import pandas as pd
from matplotlib.lines import Line2D
import numpy as np

ModelDir.set_base(__file__)

switch = 6
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
                 label="simulation(" + str(g) + "A)")
    plt.scatter(alpha100, T100, c = "b",  label="meauseremnt(" + "100" + "A)")
    plt.scatter(alpha150, T150, c = "orange",  label="meauseremnt(" + "150" + "A)")
    plt.scatter(alpha200, T200, c = "g",  label="meauseremnt(" + "200" + "A)")
    plt.scatter(alpha250, T250, c = "r",  label="meauseremnt(" + "250" + "A)")
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
    res = pd.read_pickle(ModelDir.DATA / "df_avg.pkl")

    fig = plt.figure(figsize=(6, 4))

    plt.plot(res["rotorangle"].iloc[0:451], res["tavg"].iloc[0:451], label="avg(i=100A)", c="b")
    plt.plot(res["rotorangle"].iloc[0:451], res["tmax"].iloc[0:451], label="max(i=100A)", c="b", linestyle="--")
    plt.plot(res["rotorangle"].iloc[0:451], res["tmin"].iloc[0:451], label="min(i=100A)", c="b", linestyle=(0, (1, 1)))
    plt.plot(res["rotorangle"].iloc[451:902], res["tavg"].iloc[451:902], label="avg(i=150A)", c="g")
    plt.plot(res["rotorangle"].iloc[451:902], res["tmax"].iloc[451:902], label="max(i=150A)", c="g",
             linestyle="--")
    plt.plot(res["rotorangle"].iloc[451:902], res["tmin"].iloc[451:902], label="min(i=150A)", c="g",
             linestyle=(0, (1, 1)))
    #plt.plot(res["rotorangle"].iloc[902:1353], res["tavg"].iloc[902:1353], label="avg(i=200A)", c="g")
    #plt.plot(res["rotorangle"].iloc[902:1353], res["tmax"].iloc[902:1353], label="max(i=200A)", c="g", linestyle="--")
    #plt.plot(res["rotorangle"].iloc[902:1353], res["tmin"].iloc[902:1353], label="min(i=200A)", c="g", linestyle=(0, (1, 1)))
    plt.plot(res["rotorangle"].iloc[1353:1804], res["tavg"].iloc[1353:1804], label="avg(i=250A)", c="r")
    plt.plot(res["rotorangle"].iloc[1353:1804], res["tmax"].iloc[1353:1804], label="max(i=250A)", c="r", linestyle="--")
    plt.plot(res["rotorangle"].iloc[1353:1804], res["tmin"].iloc[1353:1804], label="min(i=250A)", c="r", linestyle=(0, (1, 1)))
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


if switch == 7:
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
    res = pd.read_pickle(ModelDir.DATA / "df_avg.pkl")

    fig = plt.figure(figsize=(6, 4))
    ax1 = fig.add_subplot()
    ax1.plot(res["rotorangle"].iloc[0:451], res["tavg"].iloc[0:451], label= "i=100A")
    ax1.plot(res["rotorangle"].iloc[451:902], res["tavg"].iloc[451:902], label= "i=150A")
    ax1.plot(res["rotorangle"].iloc[902:1353], res["tavg"].iloc[902:1353], label= "i=200A")
    ax1.plot(res["rotorangle"].iloc[1353:1804], res["tavg"].iloc[1353:1804], label= "i=250A")
    plt.xlabel('Electrical angle [deg]', fontsize=10)
    plt.ylabel('Torque [Nm]', fontsize=10)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(0, 200, step=20))
    plt.yticks(fontsize=10)
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "PEMC_T7.png", bbox_inches="tight", dpi=650)
    plt.show()
