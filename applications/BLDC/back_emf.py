from matplotlib.pyplot import legend
from adze_modeler.utils import csv_read
from model import BLDCMotor

class BackEMF(BLDCMotor):
    """"""
    def __init__(self, *args, **kwargs):
        super(BackEMF, self).__init__(*args, **kwargs)

        self.msh_size_air = 1
        self.msh_size_coils = 1
        self.msh_size_magnets = 1
        self.msh_size_rotor_steel = 1
        self.msh_size_stator_steel = 1

    def add_postprocessing(self):
        super(BackEMF, self).add_postprocessing()
        # points = [(0,35.4)] # B
        points = [(-10,35.4)] # C
        # points = [(10,35.4)] # A

        self.snapshot.add_postprocessing("integration", points, "Flux")


def execute_model(m:BackEMF):
    """TODO: Docstring for execute_model.

    Args:
        m (BackEMF): TODO

    Returns: TODO

    """
    res = m(cleanup=True)
    return res['Flux']*8

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    from multiprocessing import Pool
    from adze_modeler.modelpaths import ModelDir
    from adze_modeler.utils import *

    setup_matplotlib()
    ModelDir.set_base(__file__)

    # theta = np.linspace(0, 90, 101)
    # models = [BackEMF(rotorangle=ti) for ti in theta]
    # with Pool(processes=4) as pool:
    #     res = pool.map(execute_model, models)
    #     csv_write(ModelDir.DATA/'backemf_c.csv', ['theta', 'flux'], theta, res)

    theta, flux_a = csv_read(ModelDir.DATA/'backemf_a.csv')
    theta, flux_b = csv_read(ModelDir.DATA/'backemf_b.csv')
    theta, flux_c = csv_read(ModelDir.DATA/'backemf_c.csv')

    flux_ab = np.array(flux_a) - np.array(flux_b)
    flux_ac = np.array(flux_a) - np.array(flux_c)
    flux_bc = np.array(flux_b) - np.array(flux_c)

    # plt.plot(theta, flux_a, 'r-', label='A')
    # plt.plot(theta, flux_b, 'g-', label='B')
    # plt.plot(theta, flux_c, 'b-', label='C')
    # plt.legend()
    # plt.show()


    t, Vref = csv_read(ModelDir.DATA/'backemf_ref.csv')
    t = np.array(t)
    plt.figure()
    plt.plot(np.interp(theta, [0, 90], [-0.000005476451259583872, 0.0073001095290251925]), 8*np.gradient(flux_ac, theta)/5e-6*4500, 'b-', label='simulation')
    plt.plot(t-0.0001, Vref, 'r-', label='reference')
    plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.xlabel('Time [s]')
    plt.ylabel('e.m.f [V]')
    plt.legend()
    plt.savefig(ModelDir.MEDIA / "backemf.pdf", bbox_inches="tight")
    plt.show()
