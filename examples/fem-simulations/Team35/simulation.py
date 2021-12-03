
from model import DistributedWinding

from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.encapsulator import Encapsulator, mounts
from digital_twin_distiller.simulationproject import sim
from operator import itemgetter
from math import inf


def execute_model(model: DistributedWinding):
    try:
        result = model(timeout=30, cleanup=True)
        Bz = map(itemgetter(2), result["Bz"])
        return Bz
    except:
        return None


@sim.register('default')
def default_simulation(model, modelparams, simparams, miscparams):
    x = simparams['x']
    B0 = simparams['B0']
    model = DistributedWinding(x)
    Bz = execute_model(model)
    if Bz is not None:
        f1 = max(map(lambda Bz_i: abs(Bz_i - B0), Bz))
        return {'f1': f1}
    else:
        return {'f1': inf}

if __name__ == "__main__":

    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(DistributedWinding)

    model = Encapsulator(sim)
    model.build_docs()
    mounts()
    model.run()
