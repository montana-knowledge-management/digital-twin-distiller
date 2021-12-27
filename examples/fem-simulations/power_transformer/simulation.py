from multiprocessing import Pool

from model import power_transformer

from digital_twin_distiller.encapsulator import Encapsulator
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import sim


def execute_model(model: power_transformer):
    result = model(timeout=2000, cleanup=True)
    return result


@sim.register('default')
def default_simulation(model, modelparams, simparams, miscparams):
    return "Hello World!"

if __name__ == "__main__":

    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(power_transformer)

    model = Encapsulator(sim)
    model.build_docs()
    #mounts()
    model.run()
