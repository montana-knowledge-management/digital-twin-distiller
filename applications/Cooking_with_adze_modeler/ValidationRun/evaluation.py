from matplotlib import pyplot as plt

from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from problem_solver import build
from adze_modeler.metadata import Agros2DMetadata, FemmMetadata
from random import seed, uniform
from numpy import linspace, zeros
from pathlib import Path

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

def evaluate(X):
    femm_metadata = FemmMetadata()
    femm_metadata.problem_type = "magnetic"
    femm_metadata.coordinate_type = "axisymmetric"
    femm_metadata.file_script_name = "femm_solver_script"
    femm_metadata.file_metrics_name = "femm_solution.csv"
    femm_metadata.unit = "millimeters"
    femm_metadata.smartmesh = False
    platform_femm = Femm(femm_metadata)

    agros_metadata = Agros2DMetadata()
    agros_metadata.file_script_name = "agros_solver_script"
    agros_metadata.file_metrics_name = "agros_solution.csv"
    agros_metadata.problem_type = "magnetic"
    agros_metadata.coordinate_type = "axisymmetric"
    agros_metadata.analysis_type = "steadystate"
    agros_metadata.unit = 1e-3
    agros_metadata.nb_refinements = 0
    agros_metadata.polyorder = 2
    agros_metadata.adaptivity = "hp-adaptivity"
    agros_metadata.adaptivity_tol = 1
    agros_metadata.adaptivity_steps = 100
    platform_agros = Agros2D(agros_metadata)

    platform = platform_femm

    res = build(platform, X, customid='bestof1', cleanup=False)

    Bz = [pointvalue[2] for pointvalue in res['Bz']]  # [x, y, Bz(x, y)]
    Br = [pointvalue[2] for pointvalue in res['Br']]  # [x, y, Br(x, y)]
    x = [pointvalue[0] for pointvalue in res['Br']]  # [x, y, Br(x, y)]
    y = [pointvalue[1] for pointvalue in res['Br']]  # [x, y, Br(x, y)]

    B0 = 2.0e-3
    F1 = max(map(lambda Bz_i: abs(Bz_i - B0), Bz))
    print(F1)

    returnmatrix = []

    for xi, yi, Bri, Bzi in zip(x, y, Br, Bz):
        returnmatrix.append((xi, yi, Bri, Bzi))

    return returnmatrix

def get_line_data(X):
    idxX = 0
    idxY = 1
    idxBr = 2
    idxBz = 3

    M = evaluate(X)
    M.sort(key=lambda row: row[idxX])
    scale = 1
    # scale = 1e-3


    Mline = [row for row in M if row[idxY] > 4.9*scale]


    xp = [row[idxX]*1000 for row in Mline]
    yBrp = [row[idxBr]*1000 for row in Mline]
    yBzp = [row[idxBz]*1000 for row in Mline]

    Mline = [row for row in M if (row[idxY] < 2.501*scale) and (row[idxY] > 2.49*scale)]

    xm = [row[idxX] * 1000 for row in Mline]
    yBrm = [row[idxBr] * 1000 for row in Mline]
    yBzm = [row[idxBz] * 1000 for row in Mline]

    Mline = [row for row in M if (row[idxY] < 0.01*scale) ]

    x0 = [row[idxX] * 1000 for row in Mline]
    yBr0 = [row[idxBr] * 1000 for row in Mline]
    yBz0 = [row[idxBz] * 1000 for row in Mline]

    if scale > 0.9:
        xp = [xpi / 1000 for xpi in xp]
        xm = [xmi / 1000 for xmi in xm]
        x0 = [x0i / 1000 for x0i in x0]

    return (xp, yBrp, yBzp), (x0, yBr0, yBz0), (xm, yBrm, yBzm)

if __name__ == '__main__':
    seed(42)
    X = [uniform(5.5, 50) for _ in range(4)]
    X.extend([uniform(1, 50) for _ in range(6)])

    X_oldf1 = [12.59463, 7.21559, 8.65518, 10.31800, 6.00144, 13.05700, 8.37318, 10.49246, 8.57676,11.15772]
    X_oldf1 = list(reversed(X_oldf1))

    X_newf1 = [11.82008, 12.17057, 15.43012, 12.38441, 4.82476, 12.61833, 7.24051, 14.57516, 9.05330,10.14532]
    X_newf1 = list(reversed(X_newf1))

    # from the dibarbara paper

    # 6.72152617299999e-05
    # X_1 = [8.08, 14.9, 6.74, 16.7, 5.45, 10.6, 11.7, 11.1, 13.6, 7.99]

    # 9.927623309000013e-05
    # X_1 = [7.71, 9.93, 16.9, 6.19, 18.3, 8.41, 7.19, 14.1, 8.91, 10.32]

    # 7.91522764400001e-05
    X_1 = [11, 8, 14, 10, 6, 11, 8, 11, 13, 13]

    # X_1 = list(reversed(X_1))


    upper, middle, lower = get_line_data(X_oldf1)

    refx = linspace(0, 5, 5)
    refBz = zeros(refx.shape[0]) + 2
    refBr = zeros(refx.shape[0])

    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(5, 7))

    ax[0].plot(refx, refBz, 'k-', label='Reference')
    ax[0].plot(upper[0], upper[2], 'r-o', label='z = 5 mm')
    ax[0].plot(middle[0], middle[2], 'b-o', label='z = 2.5 mm')
    ax[0].plot(lower[0], lower[2], 'g-o', label='z = 0 mm')
    # ax[0].grid(b=True, which='major', color='#666666', linestyle='-')
    ax[0].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
    ax[0].minorticks_on()
    ax[0].set_ylabel(r'B$_z$ [mT]')
    ax[0].legend()

    ax[1].plot(refx, refBr, 'k-', label='Reference')
    ax[1].plot(upper[0], upper[1], 'r-o', label='z = 5 mm')
    ax[1].plot(middle[0], middle[1], 'b-o', label='z = 2.5 mm')
    ax[1].plot(lower[0], lower[1], 'g-o', label='z = 0 mm')
    # ax[1].grid(b=True, which='major', color='#666666', linestyle='-')
    ax[1].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
    ax[1].minorticks_on()
    ax[1].set_ylabel(r'B$_r$ [mT]')
    ax[1].legend()


    plt.xlabel('r [mm]')
    plt.grid()
    # plt.savefig(Path(__file__).parent / 'media' / 'Bz_bestF1.png', dpi=550, bbox_inches='tight')
    plt.show()
