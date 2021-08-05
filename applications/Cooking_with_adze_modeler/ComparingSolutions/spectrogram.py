import operator
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from numpy import array
from numpy import meshgrid
from numpy import unique

path_base = Path(__file__).parent.parent
path_export_base = Path(__file__).parent / "media"
path_symmetric_data = path_base / "Symmetric" / "pareto_front_nsga2.csv"
path_asymmetric_data = path_base / "Asymmetric" / "pareto_front_nsga2.csv"


def get_line_from_file(filename):
    with open(filename) as f:
        yield from f


def get_processed_line_asym(filename):
    for line_i in get_line_from_file(filename):
        # for statistics_nsga2.csv
        # platform, *line_i = line_i.strip().split(',')
        # f1, f2, f3, nodes, *r = (float(ri) for ri in line_i)

        # for pareto front.csv
        line_i = line_i.strip().split(",")
        f1, f2, f3, *r = (float(ri) for ri in line_i)
        yield (f1, f2, f3, *r)


def get_processed_line_sym(filename):
    for line_i in get_line_from_file(filename):
        # for statistics_nsga2.csv
        # platform, *line_i = line_i.strip().split(',')
        # f1, f2, f3, nodes, *r = (float(ri) for ri in line_i)

        # for pareto_front_nsga2.csv
        line_i = line_i.strip().split(",")
        f1, f2, f3, *r = (float(ri) for ri in line_i)

        r = list(reversed(r))
        r.extend(reversed(r))
        yield (f1, f2, f3, *r)


###################################################################################################################
# Data acquisition

sym_data_generator = get_processed_line_sym(path_symmetric_data)
asym_data_generator = get_processed_line_asym(path_asymmetric_data)

data_sym = array([record for record in sym_data_generator])
data_asym = array([record for record in asym_data_generator])


print("len data_sym:", data_sym.shape[0])
print("len data_asym:", data_asym.shape[0])
print("-" * 30)

# Filtering out duplicate elements
data_asym = unique(data_asym, axis=0)
data_sym = unique(data_sym, axis=0)
print("len data_sym:", len(data_sym))
print("len data_asym:", len(data_asym))
print("-" * 30)

idxF1 = 0
idxF2 = 1
idxF3 = 2

###################################################################################################################

# Matplotlib setup
mm2inch = lambda x: 0.03937007874 * x
plt.rcParams["figure.figsize"] = mm2inch(160), mm2inch(100)
plt.rcParams["lines.linewidth"] = 1
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 14

plt.rc("font", size=BIGGER_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=SMALL_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=SMALL_SIZE)  # fontsize of the figure title

plt.style.use(["default", "seaborn-bright"])

X = data_sym[:, idxF1]
Y = data_sym[:, idxF2]
Z = data_sym[:, idxF3] * 2

fig, ax = plt.subplots(ncols=3, nrows=1, figsize=(20, 6), subplot_kw={"projection": "3d"})
ax[0].plot_trisurf(X, Y, Z, cmap=plt.cm.viridis)
ax[0].ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax[0].view_init(31, -8)
ax[0].set_xlabel("F1", fontsize=10)
ax[0].set_ylabel("F2", fontsize=10)
ax[0].set_zlabel("F3", fontsize=10)
ax[0].tick_params(axis="both", labelsize=10)
ax[0].set_title("a)", y=-0.2)

ax[1].plot_trisurf(X, Y, Z, cmap=plt.cm.viridis)
ax[1].ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
ax[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax[1].view_init(58, 112)
ax[1].set_xlabel("F1", fontsize=10)
ax[1].set_ylabel("F2", fontsize=10)
ax[1].set_zlabel("F3", fontsize=10)
ax[1].tick_params(axis="both", labelsize=10)
ax[1].set_title("b)", y=-0.2, fontsize=12)

ax[2].plot_trisurf(X, Y, Z, cmap=plt.cm.viridis)
ax[2].ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
ax[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax[2].view_init(23, -136)
ax[2].set_xlabel("F1", fontsize=10)
ax[2].set_ylabel("F2", fontsize=10)
ax[2].set_zlabel("F3", fontsize=10)
ax[2].tick_params(axis="both", labelsize=10)
ax[2].set_title("c)", y=-0.2)
fig.set_tight_layout(True)
plt.savefig("3in1.png", dpi=550, bbox_inches="tight")
# plt.show()
