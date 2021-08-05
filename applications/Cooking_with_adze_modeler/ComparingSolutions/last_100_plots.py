import operator
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from numpy import array
from numpy import unique

path_base = Path(__file__).parent.parent
path_export_base = Path(__file__).parent / "media"
path_520 = path_base / "Symmetric" / "pareto_front_520.csv"
path_unbounded = path_base / "Symmetric" / "pareto_front_nsga2.csv"


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

data_520_generator = get_processed_line_sym(path_520)
data_unbounded_generator = get_processed_line_sym(path_unbounded)
data_reduced_generator = get_processed_line_sym(path_reduced)

data_520 = array([record for record in data_520_generator])
data_unbounded = array([record for record in data_unbounded_generator])
data_reduced = array([record for record in data_reduced_generator])

print("len data_520:", data_520.shape[0])
print("len data_asym:", data_unbounded.shape[0])
print("len data_reduced:", data_reduced.shape[0])
print("-" * 30)

# Filtering out duplicate elements
data_unbounded = unique(data_unbounded, axis=0)
data_520 = unique(data_520, axis=0)
data_reduced = unique(data_reduced, axis=0)
print("len data_sym:", len(data_520))
print("len data_asym:", len(data_unbounded))
print("len data_reduced:", len(data_reduced))
print("-" * 30)

data_520 = data_520[-100:]
data_unbounded = data_unbounded[-100:]
data_reduced = data_reduced[-100:]

idxF1 = 0
idxF2 = 1
idxF3 = 2

###################################################################################################################

# Matplotlib setup
mm2inch = lambda x: 0.03937007874 * x
plt.rcParams["figure.figsize"] = mm2inch(40), mm2inch(40)
plt.rcParams["lines.linewidth"] = 1
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 14

plt.rc("font", size=BIGGER_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=BIGGER_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=BIGGER_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=BIGGER_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.style.use(["default", "seaborn-bright"])

###################################################################################################################
# Plot Symmetric fitness functions
data_520[:, idxF3] = data_520[:, idxF3] * 2.0
data_unbounded[:, idxF3] = data_unbounded[:, idxF3] * 2.0
data_reduced[:, idxF3] = data_reduced[:, idxF3] * 2.0

constrained = min(data_520, key=operator.itemgetter(idxF1))[idxF1]
unbounded = min(data_unbounded, key=operator.itemgetter(idxF1))[idxF1]
reduced = min(data_reduced, key=operator.itemgetter(idxF1))[idxF1]
print(f"{reduced-unbounded=:.3e}")
print(f"{reduced-constrained=:.3e}")
print(f"{reduced-unbounded=:.3e}")
print("--- " * 10)

constrained = min(data_520, key=operator.itemgetter(idxF2))[idxF2]
unbounded = min(data_unbounded, key=operator.itemgetter(idxF2))[idxF2]
reduced = min(data_reduced, key=operator.itemgetter(idxF2))[idxF2]
print(f"{reduced-unbounded=:.3e}")
print(f"{reduced-constrained=:.3e}")
print(f"{reduced-unbounded=:.3e}")
print("--- " * 10)

constrained = min(data_520, key=operator.itemgetter(idxF3))[idxF3]
unbounded = min(data_unbounded, key=operator.itemgetter(idxF3))[idxF3]
reduced = min(data_reduced, key=operator.itemgetter(idxF3))[idxF3]
print(f"{reduced-unbounded=:.3f}")
print(f"{reduced-constrained=:.3f}")
print(f"{reduced-unbounded=:.3f}")
print("--- " * 10)

data_520 = array(sorted(data_520, key=operator.itemgetter(idxF3)))
data_unbounded = array(sorted(data_unbounded, key=operator.itemgetter(idxF3)))
data_reduced = array(sorted(data_reduced, key=operator.itemgetter(idxF3)))
plt.figure()
plt.plot(data_520[:, idxF3], "r")
plt.plot(data_unbounded[:, idxF3], "g")
plt.plot(data_reduced[:, idxF3], "b")
plt.show()

# # F1 - F2
# plt.figure()
# plt.scatter(data_520[:, idxF1], data_520[:, idxF2], c='b', label='Constrained')
# plt.scatter(data_unbounded[:, idxF1], data_unbounded[:, idxF2], c='r', label='Unbounded')
# plt.scatter(data_reduced[:, idxF1], data_reduced[:, idxF2], c='g', label='Reduced')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# plt.xlabel(r"F$_1$")
# plt.ylabel(r"F$_2$")
# plt.grid(b=True, which='major', color='#666666', linestyle='-')
# plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.4)
# plt.minorticks_on()
# plt.legend()
# plt.savefig(path_export_base / 'boundedvsunbounded_last100-f1-f2_reduced.png', dpi=550, bbox_inches='tight')
# # plt.show()
# #
# # # F1 - F3
# plt.figure()
# plt.scatter(data_520[:, idxF1], data_520[:, idxF3], c='b', label='Constrained')
# plt.scatter(data_unbounded[:, idxF1], data_unbounded[:, idxF3], c='r', label='Ubounded')
# plt.scatter(data_reduced[:, idxF1], data_reduced[:, idxF3], c='g', label='Reduced')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.xlabel(r"F$_1$")
# plt.ylabel(r"F$_3$")
# plt.grid(b=True, which='major', color='#666666', linestyle='-')
# plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.4)
# plt.minorticks_on()
# plt.legend()
# plt.savefig(path_export_base / 'boundedvsunbounded_last100-f1-f3_reduced.png', dpi=550, bbox_inches='tight')
# # plt.show()
# #
# # F2 - F3
# plt.figure()
# plt.scatter(data_520[:, idxF2], data_520[:, idxF3], c='b', label='Constrained')
# plt.scatter(data_unbounded[:, idxF2], data_unbounded[:, idxF3], c='r', label='Ubounded')
# plt.scatter(data_reduced[:, idxF2], data_reduced[:, idxF3], c='g', label='Reduced')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.xlabel(r"F$_2$")
# plt.ylabel(r"F$_3$")
# plt.grid(b=True, which='major', color='#666666', linestyle='-')
# plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.4)
# plt.minorticks_on()
# plt.legend()
# plt.savefig(path_export_base / 'boundedvsunbounded_last100-f2-f3_reduced.png', dpi=550, bbox_inches='tight')
# # plt.show()
#
### Violinplot
df = {f"R{i+1}": list() for i in range(20)}
df["type"] = list()

for i in range(3, 23):
    ri = list(data_520[:, i])
    df[f"R{i+1-3}"] = ri.copy()
df["type"].extend(["Constrained"] * len(data_520))

for i in range(3, 23):
    ri = list(data_unbounded[:, i])
    df[f"R{i+1-3}"].extend(ri.copy())
df["type"].extend(["Unbounded"] * len(data_unbounded))


df = pd.DataFrame(df)
df = df.melt(value_vars=[f"R{i+1}" for i in range(20)], id_vars="type", var_name="radius_name", value_name="value")

fig, ax = plt.subplots(figsize=(5, 14))
ax = sns.violinplot(
    x="value",
    y="radius_name",
    hue="type",
    data=df,
    split=True,
    palette=["lightblue", "#F08080"],
    cut=0,
    inner=None,
    linewidth=0.5,
    width=1,
)

ax.set_xlim(0, 20)
ax.set_xlabel("Radius [mm]")
ax.set_ylabel("ith Coil")
ax.grid(b=True, which="major", color="#666666", linestyle="-")
ax.grid(b=True, which="minor", color="#999999", linestyle="-", alpha=0.4)
ax.minorticks_on()
plt.savefig(path_export_base / "boundedvsunbounded_last100-radius-distribution.png", dpi=550, bbox_inches="tight")
plt.show()


from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(16, 16))
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)
g = ax.scatter(
    data_520[:, idxF1],
    data_520[:, idxF2],
    data_520[:, idxF3] * 2,
    c="b",
    marker="o",
    depthshade=True,
    label="Constrained",
)
g = ax.scatter(
    data_unbounded[:, idxF1],
    data_unbounded[:, idxF2],
    data_unbounded[:, idxF3],
    c="r",
    marker="o",
    depthshade=True,
    label="Unbounded",
)
g = ax.scatter(
    data_reduced[:, idxF1],
    data_reduced[:, idxF2],
    data_reduced[:, idxF3],
    c="g",
    marker="o",
    depthshade=True,
    label="Reduced",
)
ax.set_xlabel(r"F$_1$")
ax.set_ylabel(r"F$_2$")
ax.set_zlabel(r"F$_3$")
plt.legend()
plt.savefig(path_export_base / "last100_3d.png", dpi=330, bbox_inches="tight")
plt.show()
