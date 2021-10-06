import csv
import multiprocessing
from operator import itemgetter
from pathlib import Path
from statistics import fmean
from time import perf_counter
from uuid import uuid4

import matplotlib.pyplot as plt
from numpy import linspace

from adze_modeler.utils import csv_read, csv_write, get_polyfit, rms
from adze_modeler.utils import setup_matplotlib, get_width_height
from model import *

DIR_SAVE = DIR_DATA / 'mesh_selectivity'

class BLDCMeshSelectivity(BLDCMotor):

    def __init__(self, rotorangle: float, meshsize:float, exportname: str = None):
        super().__init__(rotorangle=rotorangle, exportname=exportname)
        self.msh_size_airgap = meshsize
        self.msh_size_rotor_steel = meshsize
        self.msh_size_magnets = meshsize

    def add_postprocessing(self):
        super(BLDCMotor, self).add_postprocessing()
        self.snapshot.add_postprocessing("mesh_info", None, None)


def get_id():
    names = DIR_DATA / "mesh_selectivity"
    names = set(names.rglob('*.csv'))
    # print(names)
    while True:
        name = str(hex(int(uuid4())))[-8:]
        if name not in names:
            return name

def execute_model(model: BLDCMotor):
    t0 = perf_counter()
    res = model(timeout=2000, cleanup=False)
    t1 = perf_counter()
    res["Torque"]*=8
    print(f"\t{abs(model.rotorangle):.2f} ° - {abs(model.alpha):.2f} °\t {res['Torque']:.3f} Nm \t {t1-t0:.2f} s")
    return res

def clear_all():
    if DIR_SAVE.exists():
        for fi in DIR_SAVE.iterdir():
            fi.unlink()

        DIR_SAVE.rmdir()

    DIR_SAVE.mkdir()

def analyze_selectivity():
    # clear_all()
    fm = open(DIR_SAVE/'meta.csv', 'w', newline='')
    meta_writer = csv.writer(fm, quoting=csv.QUOTE_NONNUMERIC)
    meta_writer.writerow(['name', 'msh_size', 'nb_elements', 'Tpp'])

    msh_sizes = linspace(0.1, 0.6, 51)
    for msh_size_i in msh_sizes:
        name = f'M-{get_id()}.csv'
        print(f'name: {name}', f'size: {msh_size_i:.4f}', sep='\t')

        theta = linspace(0, 360 / 24 / 2, 61)
        models = [BLDCMeshSelectivity(ti, msh_size_i) for ti in theta]
        with multiprocessing.Pool(processes=12) as pool:
            t0 = perf_counter()
            results = pool.map(execute_model, models)
            t1 = perf_counter()
            T = [ri['Torque'] for ri in results]
            nb_elements = int(fmean([ri['elements'] for ri in results]))
            Tpp = 2*max(T)

            csv_write(DIR_SAVE / name, ['rotorangle', 'Torque'], theta, T)
            meta_writer.writerow([name, msh_size_i, nb_elements, Tpp])
            print(f'\tTpp: {Tpp:.5f} Nm', f'Elements: {nb_elements}', f'{t1-t0:.3} s', sep='\t')

    fm.close()

def plot_results():
    setup_matplotlib()

    data = csv_read(DIR_SAVE/'meta.csv')
    data = list(zip(*data))

    # sorting results based on the number of elements
    data.sort(key=itemgetter(1), reverse=True)

    plot_cogging_torque(data)
    plot_msh_pp(data)
    plot_msh_rms(data)

def plot_cogging_torque(d):
    data_min, *d, data_max = d
    name_min, mesh_size_min, nb_element_min, Tpp_min = data_min
    names, mesh_sizes, nb_elements, Tpps = zip(*d)
    name_max, mesh_size_max, nb_element_max, Tpp_max = data_max

    fig = plt.figure()
    ax = fig.add_subplot(111)

    theta, T = csv_read(DIR_SAVE / name_min)
    x, y = get_polyfit(theta, T, N=301)
    plt.plot(x, y, 'r-', label=f'Min. ({int(nb_element_min)})', lw=1, zorder=20)

    theta, T = csv_read(DIR_SAVE / name_max)
    x, y = get_polyfit(theta, T, N=301)
    plt.plot(x, y, 'b-', label=f'Max. ({int(nb_element_max)})', lw=1, zorder=20)

    for name in names:
        theta, T = csv_read(DIR_SAVE / name)
        x, y = get_polyfit(theta, T, N=301)
        plt.plot(x, y, 'gray', alpha=0.4)

    thetasmart, Tsmart = csv_read(DIR_SAVE / "cogging_smartmesh.csv")
    x, y = get_polyfit(thetasmart, Tsmart, N=301)
    plt.plot(x, y, 'magenta', lw=1, zorder=50, label='Smartmesh')

    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Rotor angle [°]")
    plt.ylabel("Torque [Nm]")
    plt.legend(loc="upper right")
    plt.savefig(DIR_MEDIA / "msh_selectivity.pdf", bbox_inches="tight")
    plt.show()

def plot_msh_pp(d):
    names, mesh_sizes, nb_elements, Tpp = zip(*d)
    Trms = []
    # converting the number of elements into int
    nb_elements = [int(ni) for ni in nb_elements]

    # adding rms values to the data
    for name in names:
        theta, T = csv_read(DIR_SAVE / name)
        x, y = get_polyfit(theta, T, N=301)
        Trms.append(rms(y))

    color = 'green'
    alpha = 0.7

    fig, ax = plt.subplots(nrows=1, ncols=2)
    ax[0].scatter(nb_elements, Tpp, c=color, alpha=alpha, edgecolor='k', zorder=20)
    ax[0].scatter(19408,0.5883110451692927, s=45, c='red', edgecolor='k', zorder=100)
    # ax[0].plot([0, 61000], [0.655, 0.655], 'k--')
    # ax[0].plot([0, 61000], [0.61, 0.61], 'k--')
    ax[0].grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax[0].grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax[0].minorticks_on()
    ax[0].set_xlabel("Number of elements\n\na)")
    ax[0].set_ylabel("Torque peak to peak [Nm]")
    ax[0].set_xscale('log')

    ranges = [0.54, 0.6, 0.66, 0.68]
    ax[1].hist(Tpp, bins=ranges, color=color, alpha=alpha, edgecolor='k', zorder=20)
    ax[1].grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax[1].grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax[1].minorticks_on()
    ax[1].set_xticks(ranges)
    ax[1].set_xlabel("Torque peak to peak [Nm]\n\nb)")
    ax[1].set_ylabel("Number of cases")

    plt.tight_layout()
    plt.savefig(DIR_MEDIA / "msh_pp.pdf", bbox_inches="tight")
    plt.show()

def plot_msh_rms(d):
    names, mesh_sizes, nb_elements, Tpp = zip(*d)
    Trms = []
    # converting the number of elements into int
    nb_elements = [int(ni) for ni in nb_elements]

    # adding rms values to the data
    for name in names:
        theta, T = csv_read(DIR_SAVE / name)
        x, y = get_polyfit(theta, T, N=301)
        Trms.append(rms(y))

    color = 'blue'
    alpha = 0.7

    fig, ax = plt.subplots(nrows=1, ncols=2)
    ax[0].scatter(nb_elements, Trms, c=color, alpha=alpha, edgecolor='k', zorder=20)
    # ax[0].plot([0, 61000], [0.655, 0.655], 'k--')
    # ax[0].plot([0, 61000], [0.61, 0.61], 'k--')
    ax[0].grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax[0].grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax[0].minorticks_on()
    ax[0].set_xlabel("Number of elements\n\na)")
    ax[0].set_ylabel("RMS Torque [Nm]")
    ax[0].set_xscale('log')

    ranges = [0.128, 0.145, 0.156, 0.164]
    ax[1].hist(Trms, bins=ranges, color=color, alpha=alpha, edgecolor='k', zorder=20)
    ax[1].grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    ax[1].grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    ax[1].minorticks_on()
    ax[1].set_xticks(ranges)
    ax[1].set_xlabel("RMS Torque [Nm]\n\nb)")
    ax[1].set_ylabel("Number of cases")

    plt.tight_layout()
    plt.savefig(DIR_MEDIA / "msh_rms.pdf", bbox_inches="tight")
    plt.show()

if __name__ == '__main__':


    # analyze_selectivity()

    plot_results()



    # x = list(map(itemgetter(3), data))
    # plt.hist(x, bins=5)
    # plt.show()
