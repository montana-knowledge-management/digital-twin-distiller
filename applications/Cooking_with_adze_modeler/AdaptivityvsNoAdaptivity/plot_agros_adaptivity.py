import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

current_dir = Path(__file__).parent

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


F1_no_adapt = []
F2_no_adapt = []
F3_no_adapt = []
nodes_no_adapt = []

F1_adapt = []
F2_adapt = []
F3_adapt = []
nodes_adapt = []
tol_adapt = []

with open(Path(__file__).parent / "agros_no_adaptivity.csv", "r") as f:
    for line in f:
        record = line.strip().split(',')
        key = record.pop(0)
        record = [float(ri) for ri in record]

        F1_no_adapt.append(record.pop(0))
        F2_no_adapt.append(record.pop(0))
        F3_no_adapt.append(record.pop(0))
        nodes_no_adapt.append(int(record.pop(0)))


with open(Path(__file__).parent / "agros_adaptivity.csv", "r") as f:
    for line in f:
        record = line.strip().split(',')
        key = record.pop(0)
        record = [float(ri) for ri in record]

        F1_adapt.append(record.pop(0))
        F2_adapt.append(record.pop(0))
        F3_adapt.append(record.pop(0))
        nodes_adapt.append(int(record.pop(0)))
        tol_adapt.append(record.pop(0))

# nodes_adapt_f1, F1_adapt = np.unique(np.array([nodes_adapt, F1_adapt]), axis=1)
# nodes_adapt_f2, F2_adapt = np.unique(np.array([nodes_adapt, F2_adapt]), axis=1)

plt.figure()
plt.plot(nodes_no_adapt, F1_no_adapt, 'ro', label="No Adaptivity")
plt.plot(nodes_adapt, F1_adapt, 'bo', label="Adaptivity")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("number of nodes")
plt.ylabel("F1")
plt.legend()
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.savefig(current_dir / 'f1-adaptivity.png', dpi=550, bbox_inches='tight')


plt.figure()
plt.plot(nodes_no_adapt, F2_no_adapt, 'ro', label="No Adaptivity")
plt.plot(nodes_adapt, F2_adapt, 'bo', label="Adaptivity")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("number of nodes")
plt.ylabel("F2")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.legend()
plt.savefig(current_dir / 'f2-adaptivity.png', dpi=550, bbox_inches='tight')
plt.show()