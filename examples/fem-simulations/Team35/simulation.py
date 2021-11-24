from multiprocessing import Pool

from model import Team35

from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.server import Server, mounts
from digital_twin_distiller.simulationproject import sim


def execute_model(model: Team35):
    result = model(timeout=2000, cleanup=True)
    return result


@sim.register('default')
def default_simulation(model, modelparams, simparams, miscparams):
    return "Hello World!"

if __name__ == "__main__":

    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(Team35)

    model = Server(sim)
    model.build_docs()
    mounts()
    model.run()
