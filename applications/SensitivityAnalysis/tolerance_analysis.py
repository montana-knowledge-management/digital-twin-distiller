import multiprocessing
from statistics import fmean, stdev
from time import perf_counter
import h5py
from numpy import linspace
from model import *
import itertools
from types import SimpleNamespace
from adze_modeler.utils import rms, setup_matplotlib
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial as P

def analyze_tolerance():
    cases = {
        "S1": [-0.02, 0.0, 0.02],  # mm
        "S2": [-0.05, 0.0, 0.05],  # mm
        "S4": [-0.02, 0.0, 0.02],  # mm
        "R1": [-0.05, 0.0, 0.05],  # mm
        "R3": [-0.1, 0.0, 0.1],  # °
        "R4": [-0.02, 0.0, 0.02],  # mm
        "R5": [-0.02, 0.0, 0.02],  # mm
        "R6": [-0.02, 0.0, 0.02],  # mm
        "R7": [-0.02, 0.0, 0.02],  # mm
        "airgap": [-0.02, 0.0, 0.02],  # mm
    }

    case = ["R1", "S2", "R3", "R6"]
    disturbances = list(itertools.product(*[cases[ci] for ci in case]))
    for i, di in enumerate(disturbances):
        X = dict(zip(case, di))
        print(i+1, di)
        theta = linspace(360 / 48 / 2, 360/48, 51)
        models = [ParametricPriusMotor(X, rotorangle=ti) for ti in theta]

        T = []
        t0 = perf_counter()
        with multiprocessing.Pool(processes=4) as pool:
            T = pool.map(execute_model, models)
        t1 = perf_counter()
        print(t1-t0)

        # rmtree(DIR_BASE / "snapshots")

        with h5py.File(DIR_DATA / 'datastore.hdf5', 'a') as f:
            name = models[0].name

            data = f.create_dataset(f'tolerance_analysis/{name}', (len(theta), 2), dtype='f')
            data[...] = list(zip(theta, T))

            data.attrs["theta"] = [theta[0], theta[-1], len(theta)]
            for key, value in X.items():
                data.attrs[key] = value

def get_polyfit(x, y):
    assert len(x)==len(y)
    N = 1001
    x_fine = linspace(min(x), max(x), N)
    maxy_ref = max(y)
    rmsy_ref = rms(y)
    cases = []
    for order in range(2, 19):
        p = P.fit(x, y, order)
        y_fine = p(x_fine)
        maxy_fine = max(y_fine)
        rmsy_fine = rms(y_fine)
        d1 = (maxy_fine - maxy_ref) / maxy_ref * 100
        d2 = (rmsy_fine - rmsy_ref) / rmsy_ref * 100
        cases.append((d1, d2, order))

    cases.sort(key=lambda ci: abs(ci[0]))
    best_order = cases[0][-1]
    p = P.fit(x, y, best_order)
    y_best = p(x_fine)
    print(f'Best: MAX: {cases[0][0]:.3f} % RMS: {cases[0][1]:.3f} % order: {best_order}')
    return x_fine, y_best

def get_record(data):
    res = SimpleNamespace()
    res.theta = data[:, 0]
    res.T = data[:, 1]
    res.maxT = max(res.T)
    res.theta_fine, res.T_fine = get_polyfit(res.theta, res.T) 
    res.rmsT = rms(res.T_fine)

    attrs = dict(data.attrs)
    res.attrs = attrs.copy()

    return res

def get_data():
    res = []
    with h5py.File(DIR_DATA/'datastore.hdf5', 'r') as f:
        for name, data in f['tolerance_analysis'].items():
            res.append(get_record(data))

    res.sort(key=lambda ri: ri.maxT)
    return res

def get_reference_data():
    with h5py.File(DIR_DATA/'datastore.hdf5', 'r') as f:
        res = get_record(f['cogging_torque/reference/ref_fine'])
    return res


def plot_tolerance():
    res = get_data()
    ref = get_reference_data()
    avgT = []
    stdT = []
    rmsT = [ri.rmsT for ri in res]
    for i in range(len(res[0].theta_fine)):
        avgT.append(fmean([ri.T_fine[i] for ri in res]))
        stdT.append(stdev([ri.T_fine[i] for ri in res]))

    plt.plot(res[0].theta_fine, res[0].T_fine, 'r-', label='Minimum')
    plt.plot(res[-1].theta_fine, res[-1].T_fine, 'b-', label='Maximum')
    plt.plot(res[-1].theta_fine, avgT, 'k--', label='Average')
    plt.plot(ref.theta_fine, ref.T_fine, 'magenta', label='Reference')
    plt.fill_between(res[0].theta_fine, res[-1].T_fine, res[0].T_fine,
    color='gray', alpha=0.2)
    # for ri in r:
    #     plt.plot(ri.theta_fine, ri.T_fine, 'gray')
    #     plt.plot(ri.theta, ri.T, 'ro')

    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Electrical angle [°]")
    plt.ylabel("Cogging Torque [Nm]")
    plt.legend()
    plt.savefig(DIR_MEDIA / "tolerance_noodles.pdf", bbox_inches="tight")
    plt.show()

    diff = []
    print(ref.T_fine)
    for T1,T2 in zip(res[-1].T_fine, res[0].T_fine):
        diff.append(T1-T2)


    plt.figure()
    plt.plot(res[0].theta_fine, diff, label=r"$\max T - \min T$")
    plt.plot(res[0].theta_fine, stdT, label='Std. dev.')
    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("Electrical angle [°]")
    plt.ylabel("Cogging Torque [Nm]")
    plt.legend()
    plt.savefig(DIR_MEDIA / "maxmindifference.pdf", bbox_inches="tight")
    plt.show()




    plt.figure()
    plt.hist(rmsT, bins=20)
    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel("RMS Cogging Torque [Nm]")
    plt.ylabel("Probability")
    plt.savefig(DIR_MEDIA / "rmsdist.pdf", bbox_inches="tight")
    plt.show()




if __name__ == "__main__":
    setup_matplotlib()

    # analyze_tolerance()
    plot_tolerance()
