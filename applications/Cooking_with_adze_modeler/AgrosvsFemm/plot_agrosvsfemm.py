from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

current_dir = Path(__file__).parent

mm2inch = lambda x: 0.03937007874 * x
plt.rcParams["figure.figsize"] = mm2inch(160), mm2inch(100)
plt.rcParams["lines.linewidth"] = 1
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=SMALL_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=SMALL_SIZE)  # fontsize of the figure title


plt.style.use(["default", "seaborn-bright"])


Br_agros = []
Bz_agros = []
nodes_agros = []

Br_femm = []
Bz_femm = []
nodes_femm = []


with open(Path(__file__).parent / "results.csv") as f:
    for line in f:
        record = line.strip().split(",")
        key = record.pop(0)
        record = [float(ri) for ri in record]
        Br, Bz, nb_nodes = record
        if key == "agros2d":
            if nb_nodes in nodes_agros:
                continue
            else:
                Br_agros.append(Br)
                Bz_agros.append(Bz)
                nodes_agros.append(nb_nodes)

        if key == "femm":
            if nb_nodes in nodes_femm:
                continue
            else:
                Br_femm.append(Br)
                Bz_femm.append(Bz)
                nodes_femm.append(nb_nodes)


plt.figure()
plt.plot(nodes_agros, Br_agros, "bo", label="Agros2D")
plt.plot(nodes_femm, Br_femm, "ro", label="Femm")
plt.xscale("log")
# plt.yscale("log")
plt.xlabel("number of nodes")
plt.ylabel(r"B$_r$")
plt.legend()
plt.grid(b=True, which="major", color="#666666", linestyle="-")
plt.grid(b=True, which="minor", color="#999999", linestyle="-", alpha=0.2)
plt.minorticks_on()
plt.savefig(current_dir / "Br.png", dpi=550, bbox_inches="tight")
plt.show()


plt.figure()
plt.plot(nodes_agros, Bz_agros, "bo", label="Agros2D")
plt.plot(nodes_femm, Bz_femm, "ro", label="Femm")
plt.xscale("log")
# plt.yscale("log")
plt.xlabel("number of nodes")
plt.ylabel(r"B$_z$")
plt.legend()
plt.grid(b=True, which="major", color="#666666", linestyle="-")
plt.grid(b=True, which="minor", color="#999999", linestyle="-", alpha=0.2)
plt.minorticks_on()
plt.savefig(current_dir / "Bz.png", dpi=550, bbox_inches="tight")
plt.show()
