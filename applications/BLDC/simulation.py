from adze_modeler.server import Server
from adze_modeler.simulation import sim
from model import BLDCMotor


@sim.register('SIMPLE')
def simple(model, parameters):
    return {"msg": "SIMPLE works"}


if __name__=='__main__':

    # set the model for the simulation
    sim.set_model(BLDCMotor)
    model = Server(sim)
    model.run()