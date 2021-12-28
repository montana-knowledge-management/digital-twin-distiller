from model import power_transformer

from digital_twin_distiller.encapsulator import Encapsulator
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import sim

from math import pi

def execute_model(power_transformer):
    result = model(timeout=2000, cleanup=True)
    return result


@sim.register("short_circuit_impedance")
def short_circuit_impedance(model, modelparams, simparams, miscparams):
    S = 6.3e6/3
    f = 50
    Nlv = 708
    Nhv = 650
    ff = 0.85
    js = simparams["js"]
    jp = simparams["jp"]

    m = power_transformer(js=js, jp=jp)
    res = execute_model(m)
    # res = {"Energy": 256.5673046878133}
    Wm = res["Energy"]
    Alv = modelparams['w2'] * modelparams['h2']
    Ahv = modelparams['w3'] * modelparams['h3']

    Ilv = js * Alv * ff
    Ihv = jp * Ahv * ff

    Xlv = 4 * pi * f * Wm / Ilv ** 2

    res["Xpu"] = Xlv * Ilv ** 2 / (2 * S) * 2

    return res


if __name__ == "__main__":
    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(power_transformer)

    model = Encapsulator(sim)
    # model.build_docs()
    model.host = "0.0.0.0"
    model.port = 5000
    model.run()