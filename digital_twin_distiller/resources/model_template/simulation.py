from multiprocessing import Pool

from model import ${name}

from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.server import Server
from digital_twin_distiller.simulation import sim


def execute_model(model: ${name}):
    result = model(timeout=2000, cleanup=True)
    return result


@sim.register('default')
def default_simulation(model, modelparams, simparams, miscparams):
    return "Hello World!"

if __name__ == "__main__":

    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(${name})

    model = Server(sim)
    # model.build_docs()
    model.run()
