import operator
from pathlib import Path
import matplotlib.pyplot as plt
from numpy import array, unique
import seaborn as sns
import pandas as pd

###################################################################################################################

# Matplotlib setup
mm2inch = lambda x: 0.03937007874 * x
plt.rcParams['figure.figsize'] = mm2inch(20), mm2inch(20)
plt.rcParams['lines.linewidth'] = 1
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 14

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.style.use(['default', 'seaborn-bright'])

###################################################################################################################

path_base = Path(__file__).parent.parent
path_export_base = Path(__file__).parent / "media"
path_symmetric_data = path_base / "Symmetric" / "pareto_front_nsga2.csv"
path_symmetric_data_reduced = path_base / "Symmetric" / "pareto_front_reduced.csv"


def get_line_from_file(filename):
    with open(filename, "r") as f:
        yield from f

def get_processed_line_sym(filename):
    for line_i in get_line_from_file(filename):
        # for statistics_nsga2.csv
        # platform, *line_i = line_i.strip().split(',')
        # f1, f2, f3, nodes, *r = (float(ri) for ri in line_i)

        # for pareto_front_nsga2.csv
        line_i = line_i.strip().split(',')
        f1, f2, f3, *r = (float(ri) for ri in line_i)

        r = list(reversed(r))
        r.extend(reversed(r))
        yield (f1, f2, f3, *r)

###################################################################################################################
# Data acquisition

sym_data_generator = get_processed_line_sym(path_symmetric_data)
sym_data_generator_reduced = get_processed_line_sym(path_symmetric_data_reduced)


data_sym = array([record for record in sym_data_generator])
data_sym_reduced = array([record for record in sym_data_generator_reduced])


print("len data_sym:", data_sym.shape[0])
print("len data_sym_reduced:", data_sym_reduced.shape[0])
print("-"*30)

# Filtering out duplicate elements
data_sym = unique(data_sym, axis=0)
data_sym_reduced = unique(data_sym_reduced, axis=0)
print("len data_sym:", len(data_sym))
print("len data_sym_reduced:", len(data_sym_reduced))
print("-"*30)

idxF1 = 0
idxF2 = 1
idxF3 = 2

# Plot Symmetric fitness functions
# Sort the resulst based on F1

data_sym[:, idxF3] = data_sym[:, idxF3] * 2.0
data_sym_reduced[:, idxF3] = data_sym_reduced[:, idxF3] * 2.0

data_sym = array(sorted(data_sym, key=operator.itemgetter(idxF3)))
data_sym_reduced = array(sorted(data_sym_reduced, key=operator.itemgetter(idxF3)))

N = 100
best_sym = data_sym[:N]
best_reduced = data_sym_reduced[:N]


# F1 - F2
plt.figure()
plt.scatter(best_sym[:, idxF1], best_sym[:, idxF2], c='b', label='Original')
plt.scatter(best_reduced[:, idxF1], best_reduced[:, idxF2], c='r', label='Reduced')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlabel(r"F$_1$")
plt.ylabel(r"F$_2$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
plt.minorticks_on()
plt.legend()
# plt.title(f'Best {N} solutions based on F1 sort')
#plt.savefig(path_export_base / 'bestin-F1-f1-f2.svg', format="svg", bbox_inches='tight')
plt.savefig(path_export_base / 'symmetric_bestin-F3-f1-f2_reduced.png', dpi=550, bbox_inches='tight')
# plt.show()

# F1 - F3
plt.figure()
plt.scatter(best_sym[:, idxF1], best_sym[:, idxF3], c='b', label='Original')
plt.scatter(best_reduced[:, idxF1], best_reduced[:, idxF3], c='r', label='Reduced')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"F$_1$")
plt.ylabel(r"F$_3$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
plt.minorticks_on()
plt.legend()
# plt.title(f'Best {N} solutions based on F1 sort')
#plt.savefig(path_export_base / 'bestin-F1-f1-f3.svg', format="svg", bbox_inches='tight')
plt.savefig(path_export_base / 'symmetric_bestin-F3-f1-f3_reduced.png', dpi=550, bbox_inches='tight')
# plt.show()

# F2 - F3
plt.figure()
plt.scatter(best_sym[:, idxF2], best_sym[:, idxF3], c='b', label='Original')
plt.scatter(best_reduced[:, idxF2], best_reduced[:, idxF3], c='r', label='Reduced')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"F$_2$")
plt.ylabel(r"F$_3$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
plt.minorticks_on()
plt.legend()
# plt.title(f'Best {N} solutions based on F1 sort')
#plt.savefig(path_export_base / 'bestin-F1-f2-f3.svg', format="svg", bbox_inches='tight')
plt.savefig(path_export_base / 'symmetric_bestin-F3-f2-f3_reduced.png', dpi=550, bbox_inches='tight')
# plt.show()
