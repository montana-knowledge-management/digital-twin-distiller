from math import isclose
from multiprocessing import Pool
import operator as op
import string
from time import perf_counter

import matplotlib.pyplot as plt
from numpy import linspace

from adze_modeler.doe import doe_fullfact
from adze_modeler.doe import *
from adze_modeler.utils import *
from adze_modeler.modelpaths import *
from model import BLDCMotor, execute_model

ModelDir.set_base(__file__)

DIR_SAVE = ModelDir.DATA / "doe"

class DOEBLDCMotor(BLDCMotor):

    def __init__(self, rotorangle, exportname=None, **kwargs):
        super(DOEBLDCMotor, self).__init__(rotorangle=rotorangle, exportname=exportname)
        self.airgap += kwargs.get('dairgap', 0.0)
        self.r3 += kwargs.get('dmagnet_h', 0.0)
        self.mw += kwargs.get('dmagnet_w', 0.0)
        self.dHc= kwargs.get('dHc', 0.0)
        self.dmur = kwargs.get('dmur', 0.0)


        # upodating the geometry variables
        self.r4 = self.r3+(self.airgap-self.void) / 2 # Rotor + airgap slice
        self.s1 = self.r3 + self.airgap # Stator Inner Radius
        self.s2 = self.s1 + 21.75       # Stator Outer Radius

    def define_materials(self):
        super(DOEBLDCMotor, self).define_materials()
        self.snapshot.materials['magnet'].mu_r += self.dmur
        self.snapshot.materials['magnet'].coercivity += self.dHc

def doe_full_factorial():
    dXnames = ('dairgap', 'dmagnet_h', 'dmagnet_w', 'dHc', 'dmur')
    dXvalues = (0.05, 0.05, 0.05, 5000, 0.05)
    designs = list(doe_fullfact([3] * 5))
    with Pool(processes=12) as pool:
        for i, design_i in enumerate(designs):
            fname = f'D-{i:03}.csv'
            design_i = [di-1 for di in design_i]
            disturbances = map(op.mul, design_i, dXvalues)
            dX = {name_i:di for name_i, di in zip(dXnames, disturbances)}
            print(fname, f'{i+1}/{len(designs)} {(i+1)*100/len(designs):.1f} %')

            theta = linspace(0, 360/24/2, 61)
            models = [DOEBLDCMotor(ti, **dX) for ti in theta]
            t0 = perf_counter()
            T = pool.map(execute_model, models)
            t1 = perf_counter()
            
            csv_write(DIR_SAVE / fname, ['rotorangle', 'Torque'], theta, T)
            print(f'\t Calculation time: {t1-t0:.2f} s')
            print('Max torque:', max(T))
            print('RMS torque:', rms(T))

def filter_data(data, method):
    """
    ff - Full-factorial
    """
    doe_ff = list(doe_fullfact([3] * 5) - 1)
    matcharrays = lambda a,b: all(map(isclose, a, b))
    if method=='ff':
        data.sort(key=op.itemgetter(3))
        return data
    elif method == 'bb':
        doe_designs = doe_bbdesign(5, center=1)

    elif method == 'pb':
        doe_designs = doe_pbdesign(5)

    elif method=='ccf':
        doe_designs = doe_ccf(5)

    elif method=='taguchi':
        doe_designs = [[-0, -0, -0, -0, -0],
                   [-0, 0, 0, 0, 0],
                   [-0, 1, 1, 1, 1],
                   [0, -0, -0, 0, 0],
                   [0, 0, 0, 1, 1],
                   [0, 1, 1, -0, -0],
                   [1, -0, 0, -0, 1],
                   [1, 0, 1, 0, -0],
                   [1, 1, -0, 1, 0],
                   [-0, -0, 1, 1, 0],
                   [-0, 0, -0, -0, 1],
                   [-0, 1, 0, 0, -0],
                   [0, -0, 0, 1, -0],
                   [0, 0, 1, -0, 0],
                   [0, 1, -0, 0, 1],
                   [1, -0, 1, 0, 1],
                   [1, 0, -0, 1, -0],
                   [1, 1, 0, -0, 0]
                   ]

    ret = []
    for di in doe_designs:
        idx = next(i for i,v in enumerate(doe_ff) if matcharrays(di, v))
        ret.append(data[idx])

    ret.sort(key=op.itemgetter(3))
    return ret

def plot_cogging(data:list):
    w,h = get_width_height(type_='double', unit='inch')
    fig, ax = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True, figsize=(w,h))
    xref, yref = csv_read(ModelDir.DATA /'cogging_toruqe.csv')

    #### Full-Factorial
    dmin, *d, dmax = filter_data(data, 'ff')
    N = len(d) + 2
    z=10

    c_ff = 'dimgray'
    cf_ff = 'lightgray'
    c_doe= 'darkgoldenrod'
    cf_doe = 'navajowhite'

    for i in range(2):
        for j in range(2):
            ax[i][j].plot(dmin[0], dmin[1], color=c_ff, lw=1, label=f'Full-factorial ({N})', zorder=z+1)
            ax[i][j].plot(dmax[0], dmax[1], color=c_ff, lw=1, zorder=z+1)
            ax[i][j].fill_between(dmin[0], dmin[1], dmax[1], color=cf_ff, zorder=z)
            ax[i][j].plot(xref, yref, color='black', lw=1, label='Ideal', zorder=20)

    ## Box-Behnken
    dmin, *d, dmax = filter_data(data, 'bb')
    N = len(d) + 2
    z = 10
    ax[0][0].plot(dmin[0], dmin[1], color=c_doe, lw=1, label=f'Box-Behnken ({N})', zorder=z)
    ax[0][0].plot(dmax[0], dmax[1], color=c_doe, lw=1, zorder=z)
    ax[0][0].fill_between(dmin[0], dmin[1], dmax[1], color=cf_doe, zorder=z)
    # ax[0].set_title('a)', y=-0.3, fontsize=10)


    ## Box-Behnken
    dmin, *d, dmax = filter_data(data, 'pb')
    N = len(d) + 2
    z = 10
    ax[0][1].plot(dmin[0], dmin[1], color=c_doe, lw=1, label=f'Plackett-Burman ({N})', zorder=z)
    ax[0][1].plot(dmax[0], dmax[1], color=c_doe, lw=1, zorder=z)
    ax[0][1].fill_between(dmin[0], dmin[1], dmax[1], color=cf_doe, zorder=z)
    # ax[1].set_title('b)', y=-0.3, fontsize=10)

    ## CCF
    dmin, *d, dmax = filter_data(data, 'ccf')
    N = len(d) + 2
    z = 10
    ax[1][0].plot(dmin[0], dmin[1], color=c_doe, lw=1, label=f'CCF ({N})', zorder=z)
    ax[1][0].plot(dmax[0], dmax[1], color=c_doe, lw=1, zorder=z)
    ax[1][0].fill_between(dmin[0], dmin[1], dmax[1], color=cf_doe, zorder=z)
    # ax[2].set_title('b)', y=-0.3, fontsize=10)

    ## Taguchi
    dmin, *d, dmax = filter_data(data, 'taguchi')
    N = len(d) + 2
    z = 10
    ax[1][1].plot(dmin[0], dmin[1], color=c_doe, lw=1, label=f'Taguchi ({N})', zorder=z)
    ax[1][1].plot(dmax[0], dmax[1], color=c_doe, lw=1, zorder=z)
    ax[1][1].fill_between(dmin[0], dmin[1], dmax[1], color=cf_doe, zorder=z)
    # ax[2].set_title('b)', y=-0.3, fontsize=10)

    labels = string.ascii_lowercase
    for i in range(2):
        for j in range(2):
            ax[i][j].grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
            ax[i][j].grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
            ax[i][j].minorticks_on()
            ax[i][j].legend()
            ax[i][j].set_xlim(0, 360/24/2)
            ax[i][j].set_ylim(-0.039, 0.45)
            ax[i][j].set_xlabel(f"Rotor angle [°]\n\n{labels[i * 2 + j]})")

    for i in range(2):
        ax[i][0].set_ylabel("Cogging Torque [Nm]")
    plt.tight_layout()
    plt.savefig(ModelDir.MEDIA / "doe_methods.pdf", bbox_inches="tight")
    plt.show()

def plot_pp_dist(data):
    w,h = get_width_height(type_='double', unit='inch')
    fig, ax = plt.subplots(nrows=2, ncols=2, sharey=True, figsize=(w,h))

    d = filter_data(data, 'ff')
    N = len(d)
    Tpp_ref = list(map(op.itemgetter(3), d))
    for i in range(2):
        for j in range(2):
            ax[i][j].hist(Tpp_ref, bins=15, edgecolor='k', alpha=1,
                          color='lightblue', zorder=20, label=f'Full-factorial ({N})')
    
    d = filter_data(data, 'bb')
    N = len(d)
    Tpp_doe = list(map(op.itemgetter(3), d))
    ax[0][0].hist(Tpp_doe, bins=15, edgecolor='k', alpha=1, color='firebrick', zorder=20, label=f'Box-Behnken ({N})')

    d = filter_data(data, 'pb')
    N = len(d)
    Tpp_doe = list(map(op.itemgetter(3), d))
    ax[0][1].hist(Tpp_doe, bins=15, edgecolor='k', alpha=1, color='firebrick', zorder=20, label=f'Plackett-Burman ({N})')
    

    d = filter_data(data, 'ccf')
    N = len(d)
    Tpp_doe = list(map(op.itemgetter(3), d))
    ax[1][0].hist(Tpp_doe, bins=15, edgecolor='k', alpha=1, color='firebrick', zorder=20, label=f'CCF ({N})')

    d = filter_data(data, 'taguchi')
    N = len(d)
    Tpp_doe = list(map(op.itemgetter(3), d))
    ax[1][1].hist(Tpp_doe, bins=15, edgecolor='k', alpha=1, color='firebrick', zorder=20, label=f'Taguchi ({N})')

    labels = string.ascii_lowercase
    for i in range(2):
        for j in range(2):
            ax[i][j].grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
            ax[i][j].grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
            ax[i][j].minorticks_on()
            ax[i][j].legend()
            ax[i][j].set_xlabel(f"Peak cogging torque [Nm]\n\n{labels[i * 2 + j]})")

    for i in range(2):
        ax[i][0].set_ylabel("Number of designs")
    plt.tight_layout()
    plt.savefig(ModelDir.MEDIA / "doe_pp_dist.pdf", bbox_inches="tight")
    plt.show()

def plot_rms_dist(data):
    w,h = get_width_height(type_='double', unit='inch')
    fig, ax = plt.subplots(nrows=2, ncols=2, sharey=True, figsize=(w,h))

    d = filter_data(data, 'ff')
    N = len(d)
    Trms_ref = list(map(op.itemgetter(2), d))
    for i in range(2):
        for j in range(2):
            ax[i][j].hist(Trms_ref, bins=15, edgecolor='k', alpha=1,
                       color='lightblue', zorder=20, label=f'Full-factorial ({N})')
    
    d = filter_data(data, 'bb')
    N = len(d)
    Trms_doe = list(map(op.itemgetter(2), d))
    ax[0][0].hist(Trms_doe, bins=15, edgecolor='k', alpha=1, color='firebrick', zorder=20, label=f'Box-Behnken ({N})')

    d = filter_data(data, 'pb')
    N = len(d)
    Trms_doe = list(map(op.itemgetter(2), d))
    ax[0][1].hist(Trms_doe, bins=15, edgecolor='k', alpha=1, color='firebrick', zorder=20, label=f'Plackett-Burman ({N})')
    

    d = filter_data(data, 'ccf')
    N = len(d)
    Trms_doe = list(map(op.itemgetter(2), d))
    ax[1][0].hist(Trms_doe, bins=15, edgecolor='k', alpha=1, color='firebrick', zorder=20, label=f'CCF ({N})')

    d = filter_data(data, 'taguchi')
    N = len(d)
    Trms_doe = list(map(op.itemgetter(2), d))
    ax[1][1].hist(Trms_doe, bins=15, edgecolor='k', alpha=1, color='firebrick', zorder=20, label=f'Taguchi ({N})')


    labels = string.ascii_lowercase
    for i in range(2):
        for j in range(2):
            ax[i][j].grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
            ax[i][j].grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
            ax[i][j].minorticks_on()
            ax[i][j].legend()
            ax[i][j].set_xlabel(f"RMS cogging torque [Nm]\n\n{labels[i * 2 + j]})")

    for i in range(2):
        ax[i][0].set_ylabel("Number of designs")
    plt.tight_layout()
    plt.savefig(ModelDir.MEDIA / "doe_rms_dist.pdf", bbox_inches="tight")
    plt.show()

def doe_plot():
    setup_matplotlib()
    data = []
    for fi in DIR_SAVE.iterdir():
        x,y = csv_read(fi)
        x_, y_ = get_polyfit(x, y, N=301)
        data.append((x_, y_, rms(y_), max(y)*2, int(fi.stem[2:])))

    # plot_cogging(data)
    # plot_pp_dist(data)
    plot_rms_dist(data)

if __name__ == "__main__":
    # doe_full_factorial()
    doe_plot()
