from adze_modeler.modelpaths import ModelDir
from adze_modeler.server import Server
from adze_modeler.simulation import sim
from model import PowerTransformer
from math import pi

def execute_model(model: PowerTransformer):
    result = model(timeout=2000, cleanup=True)
    return result


@sim.register('short_circuit_impedance')
def short_circuit_impedance(model, modelparams, simparams, miscparams):
    f = 50
    Nlv = 708
    Nhv = 650
    ff = 0.85
    js = simparams['js']
    jp = simparams['jp']


    m = PowerTransformer(js=js, jp=jp)
    # res = execute_model(m)
    res = {'Energy': 562.5688689280591, "ALV": 41118, "AHV": 40139}
    Wm = res['Energy']
    Alv = res.pop('ALV')
    Ahv = res.pop('AHV')

    Ilv = js * Alv / Nlv
    Ihv = js * Ahv / Nhv

    Xlv =  4*f*Wm / (ff*Nlv*Ilv**2)
    Xhv =  4*f*Wm / (ff*Nhv*Ihv**2)

    res['Xlv'] = Xlv
    res['Xhv'] = Xhv

    return res

if __name__ == "__main__":
    
    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(PowerTransformer)

    model = Server(sim)
    model.run()
