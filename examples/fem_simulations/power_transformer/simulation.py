from examples.fem_simulations.power_transformer.model import PowerTransformer

from digital_twin_distiller.encapsulator import Encapsulator
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import sim

from math import pi


def execute_model(model: PowerTransformer):
    result = model(timeout=2000, cleanup=True)
    return result


def phase_current(sb, ub, con_fact):
    """
    :param sb: nominal power [MVA]
    :param ub: line voltage
    :param con_fact: connection factor --- 1 for delta --- sqrt(3) for star connection --- sqrt(2)/2. for zig-zag
    :return:
    """
    return sb * 1e3 / ub / 3. ** 0.5 / con_fact


def calculate_base_impedance(ub, sb, f, con_fact, Wm):
    """
    :param sb: nominal power [MVA]
    :param ub: line voltage
    :param con_fact: connection factor --- 1 for delta --- sqrt(3) for star connection --- sqrt(2)/2. for zig-zag
    :param f: network frequency
    :param Wm: the magnetic field energy in Joule
    :return: short circuit impeance
    """

    # calculating the base impedance and the base currents on the selected terminal
    zb = ub ** 2. / sb  # impedance
    ib = sb * 1e3 / ub / 3. ** 0.5 / con_fact

    # calculating the impedance from the volume integrals
    L = 2 * Wm / ib ** 2.

    return 2. * pi * f * L / zb * 100.


@sim.register("short_circuit_impedance")
def short_circuit_impedance(model, modelparams, simparams, miscparams):
    """
    The FEM model gives back the stored magnetic energy in the transformer window, this quantity used to calculate the
    short-circuit impedance (sci) of the power transformer.
    The sci was normed on the impedance base, which calculated from the parameters of the LV terminal.
    """

    # S = 6.3e6 / 3
    # f = 50
    # Nlv = 708
    # Nhv = 650
    # ff = 0.85
    js = simparams["js"]
    jp = simparams["jp"]

    m = PowerTransformer(js=js, jp=jp)
    res = execute_model(m)
    # res = {"Energy": 256.5673046878133}
    Wm = res["Energy"]
    # Alv = modelparams['w2'] * modelparams['h2']
    # Ahv = modelparams['w3'] * modelparams['h3']

    # Ilv = js * Alv * ff
    # Ihv = jp * Ahv * ff

    # Xlv = 4 * pi * f * Wm / Ilv ** 2

    # base impedance for the short circuit impedance calculation
    ub = 6.9  # voltage --- kV base voltage
    sb = 10.0  # nominal power  --- MVA
    zb = ub ** 2. / sb  # impedance
    f = 50
    ib = sb * 1e3 / ub / 3. ** 0.5 / 3. ** 0.5

    # -- calculating the impedance from the volume integrals --
    L = 2 * Wm / ib ** 2.
    sci = 2. * pi * f * L / zb * 100.

    res["Xpu"] = Xlv * Ilv ** 2 / (2 * S) * 2

    return res


if __name__ == "__main__":
    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(PowerTransformer)

    model = Encapsulator(sim)
    # model.build_docs()
    model.host = "0.0.0.0"
    model.port = 5000
    model.run()
