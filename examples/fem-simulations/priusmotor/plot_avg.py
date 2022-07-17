import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from digital_twin_distiller import ModelDir
import pandas as pd
import numpy as np
from scipy import interpolate

ModelDir.set_base(__file__)

switch = 4
if switch == 0:
    res = pd.read_pickle(ModelDir.DATA / "df_avgear100.pkl")
    a1 = 0
    a2 = 26
    a3 = 3
    b1 = 29
    b2 = 36
    colors = ["#B90276", '#50237F', '#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    print(colors)
    #fig = plt.subplots(figsize=(6, 4))
    for c, e in zip(range(a1,a2,a3), range(9)):
        plt.scatter([((res["rotorangle"])[c])[d] for d in range(b1, b2)], [((res["tavg"])[c])[d] for d in range(b1, b2)], label="A=" + str(0.5+c/10) + "mm", color=colors[e])
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tavg"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(116,141,1))
        yi = spline(list(range(116,141,1)))
        plt.plot(xi,yi, color=colors[e])
    plt.xlabel('Load angle [deg]', fontsize=12)
    plt.ylabel('Average torque [Nm]', fontsize=12)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks( np.arange(116, 142, step=2))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(ncol=3, fontsize=12, loc=8)

    #plt.savefig(ModelDir.MEDIA / "ESCO6.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 1:
    res = pd.read_pickle(ModelDir.DATA / "df_avgear250.pkl")
    a1 = 0
    a2 = 26
    a3 = 3
    b1 = 32
    b2 = 38
    colors = ["#B90276", '#50237F', '#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    print(colors)
    fig = plt.subplots(figsize=(6, 4))
    for c, e in zip(range(a1,a2,a3), range(9)):
        plt.scatter([((res["rotorangle"])[c])[d] for d in range(b1, b2)], [((res["tavg"])[c])[d] for d in range(b1, b2)], label="A=" + str(0.5+c/10) + "mm", color=colors[e])
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tavg"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        print(xdata)
        xi = list(range(128,149,1))
        yi = spline(list(range(128,149,1)))
        plt.plot(xi,yi, color=colors[e])
    plt.xlabel('Terhelési szög [°]', fontsize=12)
    plt.ylabel('Nyomaték [Nm]', fontsize=12)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks( np.arange(128, 150, step=2))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(ncol=3, fontsize=10, loc=8)

    plt.savefig(ModelDir.MEDIA / "mait9.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 2:
    res1 = pd.read_pickle(ModelDir.DATA / "df_avgear100.pkl")
    res2 = pd.read_pickle(ModelDir.DATA / "df_avgear150.pkl")
    res3 = pd.read_pickle(ModelDir.DATA / "df_avgear200.pkl")
    res4 = pd.read_pickle(ModelDir.DATA / "df_avgear250.pkl")
    a1 = 0
    a2 = 26
    a3 = 1
    b1 = 0
    b2 = 46
    # fig = plt.subplots(figsize=(6, 4))
    for c in range(a1, a2, a3):
        plt.plot([((res1["rotorangle"])[c])[d] for d in range(b1, b2)], [((res1["tavg"])[c])[d] for d in range(b1, b2)], color=('#50237F'))
        plt.plot([((res2["rotorangle"])[c])[d] for d in range(b1, b2)], [((res2["tavg"])[c])[d] for d in range(b1, b2)], color=('#B90276'))
        plt.plot([((res3["rotorangle"])[c])[d] for d in range(b1, b2)], [((res3["tavg"])[c])[d] for d in range(b1, b2)], color=('#008ECF'))
        plt.plot([((res4["rotorangle"])[c])[d] for d in range(b1, b2)], [((res4["tavg"])[c])[d] for d in range(b1, b2)], color=('#78BE20'))
    plt.xlabel('Load angle [deg]', fontsize=12)
    plt.ylabel('Average torque [Nm]', fontsize=12)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(0, 200, step=20))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    custom_lines = [Line2D([0], [0], color=('#78BE20'), lw=4),
                    Line2D([0], [0], color=('#008ECF'), lw=4),
                    Line2D([0], [0], color=('#B90276'), lw=4),
                    Line2D([0], [0], color=('#50237F'), lw=4)]
    plt.legend(custom_lines, ["i=250A", "i=200A","i=150A","i=100A"], fontsize=12)
    plt.savefig(ModelDir.MEDIA / "x.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 3:
    res = pd.read_pickle(ModelDir.DATA / "df_ripple100.pkl")
    a1 = 1
    a2 = 26
    a3 = 5
    b1 = 0
    b2 = 61
    colors = ["#B90276", '#50237F', '#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    # fig = plt.subplots(figsize=(6, 4))
    for c in range(a1, a2, a3):
        plt.scatter(res["earheight"], res["maxavg"], color=colors[4])
        xdata = res["earheight"]
        ydata = res["maxavg"]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(np.linspace(0.5,3,100))
        yi = spline(xi)
    plt.plot(xi, yi, color=colors[2], label="Maximum of average torque [Nm]")
    plt.xlabel('Parameter A [deg]', fontsize=12)
    plt.ylabel('Average torque [Nm]', fontsize=12)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(0.5, 3.5, step=0.5))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)

    plt.savefig(ModelDir.MEDIA / "ESCO5.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 4:
    res = pd.read_pickle(ModelDir.DATA / "df_ripple250.pkl")
    a1 = 1
    a2 = 26
    a3 = 5
    b1 = 0
    b2 = 61
    colors = ["#B90276", '#50237F', '#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    fig = plt.subplots(figsize=(6, 4))
    for c in range(a1, a2, a3):
        plt.scatter(res["earheight"], res["maxavg"], color=colors[2])
        xdata = res["earheight"]
        ydata = res["maxavg"]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(np.linspace(0.5,3,100))
        yi = spline(xi)
    plt.plot(xi, yi, color=colors[2], label= "Átlagnyomaték maximuma [Nm]")
    plt.xlabel('"A" paraméter [mm]', fontsize=12)
    plt.ylabel('Nyomaték [Nm]', fontsize=12)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(0.5, 3.5, step=0.5))
    plt.yticks(np.arange(315.5, 321, step=0.5))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=10)

    plt.savefig(ModelDir.MEDIA / "mait10.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 5:
    fig = plt.subplots(figsize=(6, 4))
    res = pd.read_pickle(ModelDir.DATA / "df_avgear250.pkl")
    a1 = 0
    a2 = 1
    a3 = 1
    b1 = 0
    b2 = 45
    colors = ['#50237F', '#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    print(colors)
    for c, e in zip(range(a1, a2, a3), range(9)):
        # plt.plot([((res["rotorangle"])[c])[d] for d in range(b1, b2)], [((res["tavg"])[c])[d] for d in range(b1, b2)], label="A=" + str(0.5+c/10) + "mm", color=colors[e])
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tavg"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e], label= "250A " + "- átlag")
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tmin"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e], linestyle='-.', label= "250A " + "- min")
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tmax"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e], linestyle='--', label= "250A " + "- max")

    res = pd.read_pickle(ModelDir.DATA / "df_avgear150.pkl")
    a1 = 0
    a2 = 1
    a3 = 1
    b1 = 0
    b2 = 45
    colors = ['#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    print(colors)
    for c, e in zip(range(a1, a2, a3), range(9)):
        # plt.plot([((res["rotorangle"])[c])[d] for d in range(b1, b2)], [((res["tavg"])[c])[d] for d in range(b1, b2)], label="A=" + str(0.5+c/10) + "mm", color=colors[e])
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tavg"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e], label= "150A " + "- átlag")
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tmin"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e], linestyle='-.', label= "150A " + "- max")
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tmax"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e], linestyle='--', label= "150A " + "- max")

    res = pd.read_pickle(ModelDir.DATA / "df_avgear100.pkl")
    a1 = 0
    a2 = 1
    a3 = 1
    b1 = 0
    b2 = 45
    colors = ["#B90276", '#50237F', '#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    for c, e in zip(range(a1, a2, a3), range(9)):
        # plt.plot([((res["rotorangle"])[c])[d] for d in range(b1, b2)], [((res["tavg"])[c])[d] for d in range(b1, b2)], label="A=" + str(0.5+c/10) + "mm", color=colors[e])
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tavg"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e], label="100A " + "- átlag")
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tmin"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e], linestyle='-.', label="100A " + "- min")
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tmax"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e], linestyle='--', label="100A " + "- max")
    plt.xlabel('Terhelési szög [°]', fontsize=12)
    plt.ylabel('Nyomaték [Nm]', fontsize=12)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks(np.arange(0, 200, step=20))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(ncol=1, fontsize=10, loc=2)

    plt.savefig(ModelDir.MEDIA / "mait8.png", bbox_inches="tight", dpi=650)
    plt.show()

elif switch == 6:
    res = pd.read_pickle(ModelDir.DATA / "df_avgear250.pkl")
    a1 = 0
    a2 = 1
    a3 = 1
    b1 = 0
    b2 = 45
    colors = ["#B90276", '#50237F', '#005691', "#008ECF", '#00A8B0', '#78BE20', "#006249", '#525F6B', '#000']
    print(colors)
    fig = plt.subplots(figsize=(6, 4))
    for c, e in zip(range(a1,a2,a3), range(9)):
        #plt.plot([((res["rotorangle"])[c])[d] for d in range(b1, b2)], [((res["tavg"])[c])[d] for d in range(b1, b2)], label="A=" + str(0.5+c/10) + "mm", color=colors[e])
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tavg"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0,180,1))
        yi = spline(list(range(0,180,1)))
        plt.plot(xi,yi, color=colors[e])
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tmin"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e])
        xdata = [((res["rotorangle"])[c])[d] for d in range(b1, b2)]
        ydata = [((res["tmax"])[c])[d] for d in range(b1, b2)]
        spline = interpolate.InterpolatedUnivariateSpline(xdata, ydata)
        xi = list(range(0, 180, 1))
        yi = spline(list(range(0, 180, 1)))
        plt.plot(xi, yi, color=colors[e])
    plt.xlabel('Load angle [deg]', fontsize=12)
    plt.ylabel('Average torque [Nm]', fontsize=12)
    plt.grid(visible=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(visible=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xticks( np.arange(116, 142, step=2))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(ncol=3, fontsize=12, loc=8)

    #plt.savefig(ModelDir.MEDIA / "ESCO6.png", bbox_inches="tight", dpi=650)
    plt.show()