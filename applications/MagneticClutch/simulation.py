from adze_modeler.modelpaths import ModelDir
from adze_modeler.server import Server
from adze_modeler.simulation import sim
from model import MagneticClutch
from multiprocessing import Pool

def execute_model(model: MagneticClutch):
    result = model(timeout=2000, cleanup=True)
    return result


@sim.register('default')
def default_simulation(model, modelparams, simparams, miscparams):
    return "Hello World!"

if __name__ == "__main__":
    
    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(MagneticClutch)

    model = Server(sim)
    # model.build_docs()
    model.run()