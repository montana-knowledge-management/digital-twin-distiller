from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.server import Server
from digital_twin_distiller.simulationproject import sim
from model import SRM
from multiprocessing import Pool

def execute_model(model: SRM):
    result = model(timeout=2000, cleanup=True)
    return result


@sim.register('default')
def default_simulation(model, modelparams, simparams, miscparams):
    return "Hello World!"

if __name__ == "__main__":
    
    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(SRM)

    model = Server(sim)
    model.run()
