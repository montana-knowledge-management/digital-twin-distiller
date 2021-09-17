from model import *
import h5py
from math import sqrt, cos, pi
import multiprocessing
from numpy import linspace, RankWarning
import matplotlib.pyplot as plt
from adze_modeler.utils import setup_matplotlib
from numpy.polynomial import Chebyshev as CH
from numpy.polynomial import Polynomial as P

def get_chebyshev_points(a, b, n):
    assert b > a
    x = []
    for k in range(1, n+1):
        xi = 0.5 * (a + b) + 0.5 * (b - a) * cos((2 * k - 1) / (2 * n) * pi)
        x.append(xi)
    return x


def get_rms(arr):
    N = len(arr)
    rms = 1/N * sum(map(lambda arr_i: arr_i**2, arr))
    return sqrt(rms)

def get_polyfit(x, y, nb_xsteps=None):
    assert len(x) == len(y)
    N = len(x)
    maxY = max(y)
    nb_xsteps = 2 * N + 1 if nb_xsteps is None else nb_xsteps
    rms_original = get_rms(y)
    x_fine = linspace(min(x), max(x), nb_xsteps)
    variants = []
    for order in range(1, 21):
        p = CH.fit(x, y, order)
        # p = P.fit(x, y, order)
        y_fine = p(x_fine)
        rms_fine = get_rms(y_fine)
        difference = (rms_fine - rms_original) / rms_original * 100
        # maxY_fine = max(y_fine)
        # difference = (maxY_fine - maxY) / maxY * 100
        variants.append([difference, order])

    variants.sort(key=lambda pi: abs(pi[0]))
    # print(variants)
    best_order = variants[0][1]
    p = CH.fit(x, y, best_order)
    y_best = p(x_fine)
    return x_fine, y_best

def calculate_points(N:int = 3):
    # theta = linspace(360 / 48 / 2, 360 / 48, N)
    theta = get_chebyshev_points(360 / 48 / 2, 360 / 48, N)
    models = [PriusMotor(rotorangle=ti) for ti in theta]
    with multiprocessing.Pool(processes=4) as pool, h5py.File(DIR_DATA / 'datastore.hdf5', "a") as f:
        T = pool.map(execute_model, models)
        if f'{N}' in f['cogging_torque/points_chebyshev']:
            f['cogging_torque/points_chebyshev'].pop(f'{N}')

        data = f.create_dataset(f'cogging_torque/points_chebyshev/{N}', (len(T), 2), dtype='f')
        data[...] = list(zip(theta, T))


def compute_pointspolyfit():
    with h5py.File(DIR_DATA / 'datastore.hdf5', "a") as f:
        for name, data in f['cogging_torque/points_chebyshev_fine'].items():
            theta = data[:, 0]
            T = data[:, 1]
            theta1, T1 = get_polyfit(theta, T, nb_xsteps=1002)

            if name in f['cogging_torque/points_chebyshev_fine']:
                f['cogging_torque/points_chebyshev_fine'].pop(name)

            data = f.create_dataset(f'cogging_torque/points_chebyshev_fine/{name}', (len(T1), 2), dtype='f')
            data[...] = list(zip(theta1, T1))

def calculate_reference():
    theta = linspace(360 / 48 / 2, 360 / 48, 501)
    models = [PriusMotor(rotorangle=ti) for ti in theta]
    with multiprocessing.Pool(processes=4) as pool, h5py.File(DIR_DATA / 'datastore.hdf5', "a") as f:
        T = pool.map(execute_model, models)
        if 'reference' in f['cogging_torque']:
            f['cogging_torque'].pop('reference')

        data = f.create_dataset('cogging_torque/reference/ref', (len(T), 2), dtype='f')
        data[...] = list(zip(theta, T))

    compute_referencepolyfit()

def compute_referencepolyfit():
    with h5py.File(DIR_DATA / 'datastore.hdf5', "a") as f:
        data = f['cogging_torque/reference/ref']
        theta = data[:, 0]
        T = data[:, 1]
        theta1, T1 = get_polyfit(theta, T)
        print(get_rms(T), get_rms(T1))

        if 'ref_fine' in f['cogging_torque/reference']:
            f['cogging_torque/reference'].pop('ref_fine')

        data = f.create_dataset('cogging_torque/reference/ref_fine', (len(T1), 2), dtype='f')
        data[...] = list(zip(theta1, T1))

def plot():
    with h5py.File(DIR_DATA / 'datastore.hdf5', "a") as f:
        data_ref = f['cogging_torque/reference/ref']
        data_fine = f['cogging_torque/reference/ref_fine']

        theta_ref = data_ref[:, 0]
        T_ref = data_ref[:, 1]
        theta_fine = data_fine[:, 0]
        T_fine = data_fine[:, 1]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        axins = ax.inset_axes([0.1, 0.55, 0.4, 0.4])
        axins.plot(theta_ref, T_ref, "r--", linewidth=1)
        axins.plot(theta_fine, T_fine, "b-", linewidth=1)
        axins.set_xlim(5.63, 5.95)
        axins.set_ylim(1.688, 1.75)
        axins.set_xticklabels([])
        axins.set_yticklabels([])

        axins1 = ax.inset_axes([0.1, 0.1, 0.4, 0.4])
        axins1.plot(theta_ref, T_ref, "r--", linewidth=1)
        axins1.plot(theta_fine, T_fine, "b-", linewidth=1)
        axins1.set_xlim(6.7, 6.856)
        axins1.set_ylim(0.6343, 0.8)
        axins1.set_xticklabels([])
        axins1.set_yticklabels([])

        ax.indicate_inset_zoom(axins, antialiased=True, edgecolor='#1D2951', hatch=r'///',  alpha=0.7)
        ax.indicate_inset_zoom(axins1, antialiased=True, edgecolor='#1D2951', hatch='///', alpha=0.7)

        plt.plot(theta_ref, T_ref, 'r--', label='Simulation')
        plt.plot(theta_fine, T_fine, 'b-', label='Polyfit')
        plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
        plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
        plt.minorticks_on()
        plt.xlabel("Rotor angle [°]")
        plt.ylabel("Torque [Nm]")
        plt.legend(loc="lower right")
        plt.savefig(DIR_MEDIA / "cogging_ref_polyfit.pdf", bbox_inches="tight")
        plt.show()

def plot_points():
    with h5py.File(DIR_DATA / 'datastore.hdf5', "a") as f:
        for name, si in f['cogging_torque/points_chebyshev_fine'].items():
            theta_ref = si[:, 0]
            T_ref = si[:, 1]

            plt.plot(theta_ref, T_ref, '-')
        plt.show()

def compare_points():
    with h5py.File(DIR_DATA / 'datastore.hdf5', "a") as f:
        theta_ref = f['cogging_torque/reference/ref_fine'][:, 0]
        T_ref = f['cogging_torque/reference/ref_fine'][:, 1]
        rms_ref = get_rms(T_ref)
        max_ref = max(T_ref)
        names = list(f['cogging_torque/points_fine'])
        names.sort(key=lambda a: int(a))
        ch_rms_diffs = []
        ch_max_diffs = []
        uniform_rms_diffs = []
        uniform_max_diffs = []
        for name_i in names:
            theta_uniform = f[f'cogging_torque/points_fine/{name_i}'][:, 0]
            T_uniform = f[f'cogging_torque/points_fine/{name_i}'][:, 1]
            rms_uniform = get_rms(T_uniform)
            max_uniform = max(T_uniform)
            rms_uniform_diff = abs((rms_uniform - rms_ref) / rms_ref * 100)
            max_uniform_diff = abs((max_uniform - max_ref) / max_ref * 100)
            uniform_rms_diffs.append(rms_uniform_diff)
            uniform_max_diffs.append(max_uniform_diff)

            theta_ch = f[f'cogging_torque/points_chebyshev_fine/{name_i}'][:, 0]
            T_ch = f[f'cogging_torque/points_chebyshev_fine/{name_i}'][:, 1]
            rms_ch = get_rms(T_ch)
            max_ch = max(T_ch)
            rms_ch_diff = abs((rms_ch - rms_ref) / rms_ref * 100)
            max_ch_diff = abs((max_ch - max_ref) / max_ref * 100)
            ch_rms_diffs.append(rms_ch_diff)
            ch_max_diffs.append(max_ch_diff)

            # plt.plot(theta_uniform, T_uniform, 'gray', alpha=0.4)
            # plt.plot(theta_ch, T_ch, 'magenta', alpha=0.4)

            print('==' * 20)
            print(name_i)
            print(f'Reference - RMS: {rms_ref:.5f} Nm - MAX: {max_ref:.5f} Nm')
            print(f'Uniform   - RMS: {rms_uniform:.5f} Nm ({rms_uniform_diff:.3f} %) - MAX: {max_uniform:.5f} Nm ({max_uniform_diff:.3f} %)')
            print(f'Chebyshev - RMS: {rms_ch:.5f} Nm ({rms_ch_diff:.3f} %) - MAX: {max_ch:.5f} Nm ({max_ch_diff:.3f} %)')
            print('=='*20)
            print()

        names = [int(ni) for ni in names]
        names_l = [ni-0.25/2 for ni in names]
        names_r = [ni+0.25/2 for ni in names]
        fig, ax = plt.subplots(2, 1)
        ax[0].bar(names_l, uniform_rms_diffs, width=0.25, color='b', label='Uniform')
        ax[0].bar(names_r, ch_rms_diffs, width=0.25, color='g', label='Chebyshev')
        ax[0].plot([names[0], names[-1]], [0.5, 0.5], 'k--', label='Treshold')
        ax[0].set_yscale('log')
        ax[0].set_xticks(names)
        ax[0].set_ylabel("RMS")
        ax[0].legend()

        ax[1].bar(names_l, uniform_max_diffs, width=0.25, color='b', label='Uniform')
        ax[1].bar(names_r, ch_max_diffs, width=0.25, color='g', label='Chebyshev')
        ax[1].plot([names[0], names[-1]], [0.5, 0.5], 'k--', label='Treshold')
        ax[1].set_yscale('log')
        ax[1].set_xticks(names)
        ax[1].set_ylabel("MAX")
        ax[1].legend()
        plt.show()

if __name__=='__main__':
    setup_matplotlib()

    # calculate_reference()
    # plot()
    # for ni in range(3, 21):
    #     calculate_points(ni)
    #     print("=="*20)
    # compute_pointspolyfit()
    # plot_points()
    compare_points()