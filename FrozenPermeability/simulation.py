from multiprocessing import Pool

from model import FrozenPermeability

from digital_twin_distiller.encapsulator import Encapsulator
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import sim


def execute_model(model: FrozenPermeability):
    result = model(timeout=2000, cleanup=True)
    return result


@sim.register('default')
def default_simulation(model, modelparams, simparams, miscparams):
    return "Hello World"

@sim.register('simple')
def example_simulation(model, modelparams, simparams, miscparams):
    motor = model(**simparams)
    return execute_model(motor)

if __name__ == "__main__":

    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(FrozenPermeability)

    model = Encapsulator(sim)
    model.build_docs(ModelDir.DOCS)
    model.run()
