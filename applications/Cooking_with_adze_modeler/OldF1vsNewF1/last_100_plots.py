import operator
from pathlib import Path
import matplotlib.pyplot as plt
from numpy import array, unique
import seaborn as sns
import pandas as pd

path_base = Path(__file__).parent.parent
path_export_base = Path(__file__).parent / "media"
path_old_data = path_base / "ArtapOptimizationSymmetric" / "pareto_front.csv"
path_new_data = path_base / "ArtapOptimizationSymmetricNewF1"/ "pareto_front.csv"

def get_line_from_file(filename):
    with open(filename, "r") as f:
        yield from f
def get_processed_line_sym(filename):
    for line_i in get_line_from_file(filename):
        # for statistics.csv
        # platform, *line_i = line_i.strip().split(',')
        # f1, f2, f3, nodes, *r = (float(ri) for ri in line_i)

        # for pareto_front.csv
        line_i = line_i.strip().split(',')
        f1, f2, f3, *r = (float(ri) for ri in line_i)

        r = list(reversed(r))
        r.extend(reversed(r))
        yield (f1, f2, f3, *r)

###################################################################################################################
# Data acquisition

data_old_generator = get_processed_line_sym(path_old_data)
data_new_generator = get_processed_line_sym(path_new_data)

data_old = array([record for record in data_old_generator])
data_new = array([record for record in data_new_generator])


print("len data_sym:", data_old.shape[0])
print("len data_asym:", data_new.shape[0])
print("-"*30)

# Filtering out duplicate elements
data_new = unique(data_new, axis=0)
data_old = unique(data_old, axis=0)
print("len data_sym:", len(data_old))
print("len data_asym:", len(data_new))
print("-"*30)

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

plt.rc('font', size=SMALL_SIZE)          # contros default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title

plt.style.use(['default', 'seaborn-bright'])

###################################################################################################################
# Plot Symmetric fitness functions
# Sort the resulst based on F3

data_old[:, idxF3] = data_old[:, idxF3] * 2.0
data_new[:, idxF3] = data_new[:, idxF3] * 2.0

# plt.figure()
# plt.plot(data_sym[:, idxF1])
# plt.plot(data_sym[:, idxF2])
# plt.plot(data_sym[:, idxF3]/100000)
# plt.show()


# F1 - F2
plt.figure()
plt.scatter(data_old[:, idxF1], data_old[:, idxF2], c='b', label='Old')
plt.scatter(data_new[:, idxF1], data_new[:, idxF2], c='r', label='New')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"F$_1$")
plt.ylabel(r"F$_2$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.legend()
plt.title('Last 100 individuals')
plt.savefig(path_export_base / 'last100-f1-f2.png', dpi=550, bbox_inches='tight')
# plt.show()

# F1 - F3
plt.figure()
plt.scatter(data_old[:, idxF1], data_old[:, idxF3], c='b', label='Old')
plt.scatter(data_new[:, idxF1], data_new[:, idxF3], c='r', label='New')
# plt.scatter(rest_asym[:, idxF1], rest_asym[:, idxF3], c='g', alpha=0.5, label='Rest Asymmetric')
# plt.scatter(rest_sym[:, idxF1], rest_sym[:, idxF3], c='magenta', alpha=0.1, label='Rest Symmetric')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"F$_1$")
plt.ylabel(r"F$_3$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.legend()
plt.title('Last 100 individuals')
plt.savefig(path_export_base / 'last100-f1-f3.png', dpi=550, bbox_inches='tight')
# plt.show()

# F2 - F3
plt.figure()
plt.scatter(data_old[:, idxF2], data_old[:, idxF3], c='b', label='Old')
plt.scatter(data_new[:, idxF2], data_new[:, idxF3], c='r', label='New')
# plt.scatter(rest_asym[:, idxF2], rest_asym[:, idxF3], c='g', alpha=0.5, label='Rest Asymmetric')
# plt.scatter(rest_sym[:, idxF2], rest_sym[:, idxF3], c='magenta', alpha=0.1, label='Rest Symmetric')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"F$_2$")
plt.ylabel(r"F$_3$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.legend()
plt.title('Last 100 individuals')
plt.savefig(path_export_base / 'last100-f2-f3.png', dpi=550, bbox_inches='tight')
# plt.show()

### Violinplot
df = {f"R{i+1}": list() for i in range(20)}
df['type'] = list()

for i in range(3, 23):
    ri = list(data_new[:, i])
    df[f"R{i+1-3}"] = ri.copy()
df["type"].extend(['New']*len(data_new))

for i in range(3, 23):
    ri = list(data_old[:, i])
    df[f"R{i+1-3}"].extend(ri.copy())
df["type"].extend(['Old']*len(data_old))

df = pd.DataFrame(df)
df = df.melt(value_vars=[f"R{i+1}" for i in range(20)], id_vars='type', var_name='radius_name', value_name='value')

fig, ax = plt.subplots(figsize=(5, 14))
ax = sns.violinplot(x="value", y="radius_name",
                    hue="type",
                    data=df,
                    split=True,
                    palette=['lightblue', '#F08080'],
                    cut=0,
                    inner=None,
                    linewidth=0.5,
                    width=1)

ax.set_xlabel("Radius [mm]")
ax.set_ylabel("ith Coil")
ax.grid(b=True, which='major', color='#666666', linestyle='-')
ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
ax.minorticks_on()
plt.savefig(path_export_base/'last100-radius-distribution.png', dpi=550, bbox_inches='tight')