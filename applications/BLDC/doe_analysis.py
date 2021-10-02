import matplotlib.pyplot as plt

from model import *
from numpy import linspace
from multiprocessing import Pool
from adze_modeler.doe import fullfact
import operator as op
from time import perf_counter
from adze_modeler.utils import csv_write, csv_read, get_polyfit, setup_matplotlib, rms

DIR_SAVE = DIR_DATA / "doe"

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
    designs = list(fullfact([3]*5))
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

def plot_cogging(d:list):
    plt.figure()
    xref, yref = csv_read(DIR_DATA/'cogging_toruqe.csv')

    d.sort(key=op.itemgetter(3))
    dmin, *d, dmax = d
    for x, y, rmsT, maxT in d:
        plt.plot(x, y, 'gray', alpha=0.5)
    
    plt.plot(dmin[0], dmin[1], 'b-', lw=2, label=f'Min. peak torque', zorder=20)
    plt.plot(dmax[0], dmax[1], 'r-', lw=2, label=f'Max. peak torque', zorder=20)
    plt.plot(xref, yref, color='magenta', lw=2, label='Ideal', zorder=20)

    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.legend()
    plt.xlim(0, 360/24/2)
    plt.ylim(-0.039, 0.45)
    plt.xlabel("Rotor angle [°]")
    plt.ylabel("Cogging Torque [Nm]")
    plt.savefig(DIR_MEDIA / "doe_noodles.pdf", bbox_inches="tight")
    plt.show()

def plot_pp_dist(d):
    Tpp = tuple(map(op.itemgetter(3), d))

    plt.hist(Tpp, bins=15, edgecolor='k', alpha=1, color='lightblue', zorder=20)
    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Peak cogging torque [Nm]")
    plt.ylabel("Number of designs")
    plt.savefig(DIR_MEDIA / "doe_pp_dist.pdf", bbox_inches="tight")
    plt.show()

def plot_rms_dist(d):
    Trms = tuple(map(op.itemgetter(2), d))

    plt.hist(Trms, bins=15, edgecolor='k', alpha=1, color='lightgreen', zorder=20, density=True)
    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("RMS cogging torque [Nm]")
    plt.ylabel("Number of designs")
    plt.savefig(DIR_MEDIA / "doe_rms_dist.pdf", bbox_inches="tight")
    plt.show()

def doe_plot():
    setup_matplotlib()
    data = []
    for fi in DIR_SAVE.iterdir():
        x,y = csv_read(fi)
        x_, y_ = get_polyfit(x, y, N=301)
        data.append((x_, y_, rms(y_), max(y)))

    # plot_cogging(data)
    # plot_pp_dist(data)
    plot_rms_dist(data)

if __name__ == "__main__":
    # doe_full_factorial()
    doe_plot()
