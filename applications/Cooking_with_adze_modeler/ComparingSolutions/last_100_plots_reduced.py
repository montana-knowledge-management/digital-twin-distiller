import operator
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from numpy import array
from numpy import unique

path_base = Path(__file__).parent.parent
path_export_base = Path(__file__).parent / "media"
path_symmetric_data = path_base / "Symmetric" / "pareto_front_reduced.csv"
path_asymmetric_data = path_base / "Asymmetric" / "pareto_front_reduced.csv"


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
# data_asym = unique(data_asym, axis=0)
# data_sym = unique(data_sym, axis=0)
print("len data_sym:", len(data_sym))
print("len data_asym:", len(data_asym))
print("-" * 30)

data_sym = data_sym[-100:]
data_asym = data_asym[-100:]

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
data_sym[:, idxF3] = data_sym[:, idxF3] * 2.0

# plt.figure()
# plt.plot(data_sym[:, idxF1], 'r')
# plt.plot(data_sym[:, idxF2], 'g')
# plt.plot(data_sym[:, idxF3]/100000, 'b')
# plt.show()
#
# # F1 - F2
# plt.figure()
# plt.scatter(data_sym[:, idxF1], data_sym[:, idxF2], c='b', label='Symmetric')
# plt.scatter(data_asym[:, idxF1], data_asym[:, idxF2], c='r', label='Asymmetric')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
# plt.xlabel(r"F$_1$")
# plt.ylabel(r"F$_2$")
# plt.grid(b=True, which='major', color='#666666', linestyle='-')
# plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.4)
# plt.minorticks_on()
# plt.legend()
# plt.title('Last 100 individuals')
# plt.savefig(path_export_base / 'last100-f1-f2_nsga2.png', dpi=550, bbox_inches='tight')
# plt.show()
# #
# # F1 - F3
# plt.figure()
# plt.scatter(data_sym[:, idxF1], data_sym[:, idxF3], c='b', label=r'$2 \times$Symmetric')
# plt.scatter(data_asym[:, idxF1], data_asym[:, idxF3], c='r', label='Asymmetric')
# # plt.scatter(rest_asym[:, idxF1], rest_asym[:, idxF3], c='g', alpha=0.5, label='Rest Asymmetric')
# # plt.scatter(rest_sym[:, idxF1], rest_sym[:, idxF3], c='magenta', alpha=0.1, label='Rest Symmetric')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.xlabel(r"F$_1$")
# plt.ylabel(r"F$_3$")
# plt.grid(b=True, which='major', color='#666666', linestyle='-')
# plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.4)
# plt.minorticks_on()
# plt.legend()
# plt.title('Last 100 individuals')
# plt.savefig(path_export_base / 'last100-f1-f3_nsga2.png', dpi=550, bbox_inches='tight')
# plt.show()
# #
# # F2 - F3
# plt.figure()
# plt.scatter(data_sym[:, idxF2], data_sym[:, idxF3], c='b', label=r'$2 \times$Symmetric')
# plt.scatter(data_asym[:, idxF2], data_asym[:, idxF3], c='r', label='Asymmetric')
# # plt.scatter(rest_asym[:, idxF2], rest_asym[:, idxF3], c='g', alpha=0.5, label='Rest Asymmetric')
# # plt.scatter(rest_sym[:, idxF2], rest_sym[:, idxF3], c='magenta', alpha=0.1, label='Rest Symmetric')
# plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
# plt.xlabel(r"F$_2$")
# plt.ylabel(r"F$_3$")
# plt.grid(b=True, which='major', color='#666666', linestyle='-')
# plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.4)
# plt.minorticks_on()
# plt.legend()
# plt.title('Last 100 individuals')
# plt.savefig(path_export_base / 'last100-f2-f3_nsga2.png', dpi=550, bbox_inches='tight')
# plt.show()

### Violinplot
df = {f"R{i+1+4}": list() for i in range(12)}
df["type"] = list()

for i in range(3, 15):
    ri = list(data_sym[:, i])
    df[f"R{i+1-3+4}"] = ri.copy()
df["type"].extend(["Symmetric"] * len(data_sym))

for i in range(3, 15):
    ri = list(data_asym[:, i])
    df[f"R{i+1-3+4}"].extend(ri.copy())
df["type"].extend(["Asymmetric"] * len(data_asym))

df["R1"] = [0] * len(df["R5"])
df["R2"] = [0] * len(df["R5"])
df["R3"] = [0] * len(df["R5"])
df["R4"] = [0] * len(df["R5"])
df["R17"] = [0] * len(df["R5"])
df["R18"] = [0] * len(df["R5"])
df["R19"] = [0] * len(df["R5"])
df["R20"] = [0] * len(df["R5"])

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
# plt.savefig(path_export_base/'last100-radius-distribution_reduced.png', dpi=550, bbox_inches='tight')
plt.show()

# def is_symmetric(X):
#     assert len(X) == 20
#     upper = X[:10]
#     lower = X[10:]
#     for ui, li in zip(upper, lower):
#         if abs(ui-li) > 0.5:
#             return False
#
#     return True
#
# symmetric_solutions = list(filter(is_symmetric, data_asym[:, idxF3+1:]))
# print(len(symmetric_solutions))
