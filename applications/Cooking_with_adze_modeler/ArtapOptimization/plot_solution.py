import matplotlib.pyplot as plt

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

agrosdata = {
    "F1": [],
    "F2": [],
    "F3": [],
    "nodes": [],
    "r": []
}

femmdata = {
    "F1": [],
    "F2": [],
    "F3": [],
    "nodes": [],
    "r": []
}

with open("statistics.csv", "r") as f:
    for line in f:
        record = line.strip().split(',')
        key = record.pop(0)
        record = [float(ri) for ri in record]
        if key=='agros2d':
            agrosdata['F1'].append(record.pop(0))
            agrosdata['F2'].append(record.pop(0))
            agrosdata['F3'].append(record.pop(0))
            agrosdata['nodes'].append(record.pop(0))
            agrosdata['r'].append(record.copy())

        if key=='femm':
            femmdata['F1'].append(record.pop(0))
            femmdata['F2'].append(record.pop(0))
            femmdata['F3'].append(record.pop(0))
            femmdata['nodes'].append(record.pop(0))
            femmdata['r'].append(record.copy())

# plt.figure()
# plt.scatter(agrosdata["F1"], agrosdata["F2"])
# plt.show()

plt.figure()
plt.scatter(agrosdata["F1"], agrosdata["F2"], c='b', label='Agros')
# plt.scatter(femmdata["F1"], femmdata["F2"], c='r', label='Femm')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"F$_1$")
plt.ylabel(r"F$_2$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()
plt.legend()
plt.savefig('f1-f2.svg', format="svg", bbox_inches='tight')
plt.savefig('f1-f2.png', dpi=550, bbox_inches='tight')
# plt.show()

plt.figure()
plt.scatter(agrosdata["F1"], agrosdata["F3"], c='b', label='Agros')
# plt.scatter(femmdata["F1"], femmdata["F3"], c='r', label='Femm')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"F$_1$")
plt.ylabel(r"F$_3$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.minorticks_on()
plt.legend()
plt.savefig('f1-f3.svg', format="svg", bbox_inches='tight')
plt.savefig('f1-f3.png', dpi=550, bbox_inches='tight')
# plt.show()

plt.figure()
plt.scatter(agrosdata["F2"], agrosdata["F3"], c='b', label='Agros')
# plt.scatter(femmdata["F2"], femmdata["F3"], c='r', label='Femm')
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel(r"F$_2$")
plt.ylabel(r"F$_3$")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
plt.minorticks_on()
plt.legend()
plt.legend()
plt.savefig('f2-f3.svg', format="svg", bbox_inches='tight')
plt.savefig('f2-f3.png', dpi=550, bbox_inches='tight')
plt.show()