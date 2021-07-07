import operator
from pathlib import Path
import matplotlib.pyplot as plt
from numpy import array, unique
import seaborn as sns
import pandas as pd

path_base = Path(__file__).parent.parent
path_export_base = Path(__file__).parent / "media"
path_export_sym = path_export_base / "symmetric"
path_export_asym = path_export_base / "asymmetric"
path_export_comp = path_export_base / "comparison"
path_symmetric_data = path_base / "ArtapOptimizationSymmetric" / "statistics.csv"
path_asymmetric_data = path_base / "ArtapOptimization" / "statistics.csv"


def get_line_from_file(filename):
    with open(filename, "r") as f:
        yield from f


def get_processed_line_asym(filename):
    for line_i in get_line_from_file(filename):
        platform, *line_i = line_i.strip().split(',')
        f1, f2, f3, nodes, *r = (float(ri) for ri in line_i)
        yield (f1, f2, f3, *r)


def get_processed_line_sym(filename):
    for line_i in get_line_from_file(filename):
        platform, *line_i = line_i.strip().split(',')
        f1, f2, f3, nodes, *r = (float(ri) for ri in line_i)
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
plt.rcParams['figure.figsize'] = mm2inch(160), mm2inch(100)
plt.rcParams['lines.linewidth'] = 1
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title

plt.style.use(['default', 'seaborn-bright'])

###################################################################################################################
# Plot Symmetric fitness functions
# Sort the resulst based on F2
data_sym = array(sorted(data_sym, key=operator.itemgetter(idxF2)))
data_asym = array(sorted(data_asym, key=operator.itemgetter(idxF2)))

N = 100
best_sym = data_sym[:N]
rest_sym = data_sym[N:]

best_asym = data_asym[:N]
rest_asym = data_sym[N:]

# F1 - F2
plt.figure()
plt.scatter(best_sym[:, idxF1], best_sym[:, idxF2], c='b', label='Symmetric')
plt.scatter(best_asym[:, idxF1], best_asym[:, idxF2], c='r', label='Asymmetric')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
plt.xlabel(r"F$_1$")
plt.ylabel(r"F$_2$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.legend()
plt.title(f'Best {N} solutions based on F2 sort')
plt.savefig(path_export_comp / 'svg' / 'bestin-F2-f1-f2.svg', format="svg", bbox_inches='tight')
plt.savefig(path_export_comp / 'png' / 'bestin-F2-f1-f2.png', dpi=550, bbox_inches='tight')
# plt.show()

# F1 - F3
plt.figure()
plt.scatter(best_sym[:, idxF1], best_sym[:, idxF3] * 2, c='b', label=r'$2 \times$Symmetric')
plt.scatter(best_asym[:, idxF1], best_asym[:, idxF3], c='r', label='Asymmetric')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
plt.xlabel(r"F$_1$")
plt.ylabel(r"F$_3$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.legend()
plt.title(f'Best {N} solutions based on F2 sort')
plt.savefig(path_export_comp / 'svg' / 'bestin-F2-f1-f3.svg', format="svg", bbox_inches='tight')
plt.savefig(path_export_comp / 'png' / 'bestin-F2-f1-f3.png', dpi=550, bbox_inches='tight')
# plt.show()

# F2 - F3
plt.figure()
plt.scatter(best_sym[:, idxF2], best_sym[:, idxF3] * 3, c='b', label=r'$2 \times$Symmetric')
plt.scatter(best_asym[:, idxF2], best_asym[:, idxF3], c='r', label='Asymmetric')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
plt.xlabel(r"F$_2$")
plt.ylabel(r"F$_3$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.legend()
plt.title(f'Best {N} solutions based on F2 sort')
plt.savefig(path_export_comp / 'svg' / 'bestin-F2-f2-f3.svg', format="svg", bbox_inches='tight')
plt.savefig(path_export_comp / 'png' / 'bestin-F2-f2-f3.png', dpi=550, bbox_inches='tight')
# plt.show()

### Violinplot
df = {f"R{i + 1}": list() for i in range(20)}
df['type'] = list()

for i in range(3, 23):
    ri = list(best_sym[:, i])
    df[f"R{i + 1 - 3}"] = ri.copy()
df["type"].extend(['Symmetric'] * len(best_sym))

for i in range(3, 23):
    ri = list(best_asym[:, i])
    df[f"R{i + 1 - 3}"].extend(ri.copy())
df["type"].extend(['Asymmetric'] * len(best_asym))

df = pd.DataFrame(df)
df = df.melt(value_vars=[f"R{i + 1}" for i in range(20)], id_vars='type', var_name='radius_name', value_name='value')

fig, ax = plt.subplots(figsize=(10, 5))
ax = sns.violinplot(x="radius_name", y="value",
                    hue="type",
                    data=df,
                    split=True,
                    palette=['lightblue', '#F08080'],
                    cut=0,
                    inner=None,
                    linewidth=0.5,
                    scale="count")
ax.set_xlabel("ith Coil")
ax.set_ylabel("Radius [mm]")
ax.grid(b=True, which='major', color='#666666', linestyle='-')
ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
ax.minorticks_on()
plt.savefig(path_export_comp / 'svg' / 'bestin-F2-radius-distribution.svg', format="svg", bbox_inches='tight')
plt.savefig(path_export_comp / 'png' / 'bestin-F2-radius-distribution.png', dpi=550, bbox_inches='tight')
# plt.show()

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(16, 16))
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)
g = ax.scatter(best_sym[:, idxF1], best_sym[:, idxF2], best_sym[:, idxF3] * 2, c='b', marker='o', depthshade=True,
               label=r'$2 \times$Symmetric')
g = ax.scatter(best_asym[:, idxF1], best_asym[:, idxF2], best_asym[:, idxF3], c='r', marker='o', depthshade=True,
               label='Asymmetric')
ax.set_xlabel(r'F$_1$')
ax.set_ylabel(r'F$_2$')
ax.set_zlabel(r'F$_3$')
plt.legend()
plt.savefig(path_export_comp / 'png' / 'bestin-F2-3D.png', dpi=330, bbox_inches='tight')
# plt.show()
