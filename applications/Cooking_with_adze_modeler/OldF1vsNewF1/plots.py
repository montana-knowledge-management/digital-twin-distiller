from pathlib import Path
import matplotlib.pyplot as plt
from numpy import array, unique
import operator
import seaborn as sns
import pandas as pd

path_base = Path(__file__).parent.parent
path_export_base = Path(__file__).parent / "media"

path_symmetric_data_old = path_base / "SymmetricOnlyF1" / "pareto_front_nsga2.csv"
path_asymmetric_data_old = path_base / "AsymmetricOnlyF1" / "pareto_front_nsga2.csv"

path_symmetric_data_new = path_base / "SymmetricOnlyF1New" / "pareto_front_nsga2.csv"
path_asymmetric_data_new = path_base / "AsymmetricOnlyF1New" / "pareto_front_nsga2.csv"

def get_line_from_file(filename):
    with open(filename, "r") as f:
        yield from f

def get_processed_line_asym(filename):
    for line_i in get_line_from_file(filename):
        line_i = line_i.strip().split(',')
        f1, *r = (float(ri) for ri in line_i)
        yield (f1, *r)

def get_processed_line_sym(filename):
    for line_i in get_line_from_file(filename):
        line_i = line_i.strip().split(',')
        f1, *r = (float(ri) for ri in line_i)

        r = list(reversed(r))
        r.extend(reversed(r))
        yield (f1, *r)


###################################################################################################################
# Data acquisition

sym_data_old_generator = get_processed_line_sym(path_symmetric_data_old)
asym_data_old_generator = get_processed_line_asym(path_asymmetric_data_old)

sym_data_new_generator = get_processed_line_sym(path_symmetric_data_new)
asym_data_new_generator = get_processed_line_asym(path_asymmetric_data_new)

data_sym_old = array([record for record in sym_data_old_generator])
data_asym_old = array([record for record in asym_data_old_generator])

data_sym_new = array([record for record in sym_data_new_generator])
data_asym_new = array([record for record in asym_data_new_generator])


print("len data_sym_old:", data_sym_old.shape[0])
print("len data_asym_old:", data_asym_old.shape[0])
print()
print("len data_sym_new:", data_sym_new.shape[0])
print("len data_asym_new:", data_asym_new.shape[0])
print("-"*30)

# Filtering out duplicate elements
data_sym_old = unique(data_sym_old, axis=0)
data_asym_old = unique(data_asym_old, axis=0)
data_sym_new = unique(data_sym_new, axis=0)
data_asym_new = unique(data_asym_new, axis=0)

print("len data_sym_old:", data_sym_old.shape[0])
print("len data_asym_old:", data_asym_old.shape[0])
print()
print("len data_sym_new:", data_sym_new.shape[0])
print("len data_asym_new:", data_asym_new.shape[0])
print("-"*30)

idxF1 = 0
# idxF2 = 1
# idxF3 = 2

###################################################################################################################

# Matplotlib setup
mm2inch = lambda x: 0.03937007874 * x
plt.rcParams['figure.figsize'] = mm2inch(160), mm2inch(100)
plt.rcParams['lines.linewidth'] = 1
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title

plt.style.use(['default', 'seaborn-bright'])

###################################################################################################################
# Sorting
data_sym_old = array(sorted(data_sym_old, key=operator.itemgetter(idxF1)))
data_asym_old = array(sorted(data_asym_old, key=operator.itemgetter(idxF1)))

data_sym_new = array(sorted(data_sym_new, key=operator.itemgetter(idxF1)))
data_asym_new = array(sorted(data_asym_new, key=operator.itemgetter(idxF1)))


N = 100
data_sym_old_best = data_sym_old[:N]
data_asym_old_best = data_asym_old[:N]

data_sym_new_best = data_sym_new[:N]
data_asym_new_best = data_asym_new[:N]

###################################################################################################################





###################################################################################################################
# PLOTS

plt.figure()
plt.semilogy(data_sym_old_best[:, idxF1], 'b-', label='Symmetric')
plt.semilogy(data_asym_old_best[:, idxF1], 'r-', label='Asymmetric')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"index")
plt.ylabel(r"F$_1$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
plt.minorticks_on()
plt.legend()
plt.title(f'Old F1')
#plt.savefig(path_export_base / 'bestin-F1-f1-f2.svg', format="svg", bbox_inches='tight')
# plt.savefig(path_export_base / 'bestin-F1-f1-f2.png', dpi=550, bbox_inches='tight')
plt.show()


plt.figure()
plt.semilogy(data_sym_new_best[:, idxF1], 'b-', label='Symmetric')
plt.semilogy(data_asym_new_best[:, idxF1], 'r-', label='Asymmetric')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"index")
plt.ylabel(r"F$_1$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
plt.minorticks_on()
plt.legend()
plt.title(f'New F1')
#plt.savefig(path_export_base / 'bestin-F1-f1-f2.svg', format="svg", bbox_inches='tight')
# plt.savefig(path_export_base / 'bestin-F1-f1-f2.png', dpi=550, bbox_inches='tight')
plt.show()


plt.figure()
plt.semilogy(data_sym_old_best[:, idxF1], 'b-', label='Old')
plt.semilogy(data_sym_new_best[:, idxF1], 'r-', label='New')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"index")
plt.ylabel(r"F$_1$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
plt.minorticks_on()
plt.legend()
plt.title(f'Symmetric Case')
#plt.savefig(path_export_base / 'bestin-F1-f1-f2.svg', format="svg", bbox_inches='tight')
# plt.savefig(path_export_base / 'bestin-F1-f1-f2.png', dpi=550, bbox_inches='tight')
plt.show()


plt.figure()
plt.semilogy(data_asym_old_best[:, idxF1], 'b-', label='Old')
plt.semilogy(data_asym_new_best[:, idxF1], 'r-', label='New')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"index")
plt.ylabel(r"F$_1$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
plt.minorticks_on()
plt.legend()
plt.title(f'Asymmetric Case')
#plt.savefig(path_export_base / 'bestin-F1-f1-f2.svg', format="svg", bbox_inches='tight')
# plt.savefig(path_export_base / 'bestin-F1-f1-f2.png', dpi=550, bbox_inches='tight')
plt.show()

print()
print(f'Best Sym Old: {data_sym_old[0, 0]:>8.5e}' , ' '.join([f'{ri:>8.5f}' for ri in data_sym_old[0, 1:]]))
print(f'Best Sym New: {data_sym_new[0, 0]:>8.5e}', ' '.join([f'{ri:>8.5f}' for ri in data_sym_new[0, 1:]]))
print('-' * 20)
print(f'Best Asym Old: {data_asym_old[0, 0]:>8.5e}', ' '.join([f'{ri:>8.5f}' for ri in data_asym_old[0, 1:]]))
print(f'Best Asym New: {data_asym_new[0, 0]:>8.5e}', ' '.join([f'{ri:>8.5f}' for ri in data_asym_new[0, 1:]]))


### Violinplot
df = {f"R{i+1}": list() for i in range(20)}
df['type'] = list()

for i in range(1, 21):
    ri = list(data_sym_old_best[:, i])
    df[f"R{i}"] = ri.copy()
df["type"].extend(['Old']*len(data_sym_old_best))

for i in range(1, 21):
    ri = list(data_sym_new_best[:, i])
    df[f"R{i}"].extend(ri.copy())
df["type"].extend(['New']*len(data_sym_new_best))

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
                    linewidth=1,)

ax.set_xlabel("Radius [mm]")
ax.set_ylabel("ith Coil")
ax.grid(b=True, which='major', color='#666666', linestyle='-')
ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
ax.minorticks_on()
plt.show()










### Violinplot
df = {f"R{i+1}": list() for i in range(20)}
df['type'] = list()

for i in range(1, 21):
    ri = list(data_asym_old_best[:, i])
    df[f"R{i}"] = ri.copy()
df["type"].extend(['Old']*len(data_asym_old_best))

for i in range(1, 21):
    ri = list(data_asym_new_best[:, i])
    df[f"R{i}"].extend(ri.copy())
df["type"].extend(['New']*len(data_asym_new_best))

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
                    linewidth=1,)

ax.set_xlabel("Radius [mm]")
ax.set_ylabel("ith Coil")
ax.grid(b=True, which='major', color='#666666', linestyle='-')
ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
ax.minorticks_on()
plt.show()