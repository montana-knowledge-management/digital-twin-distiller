from adze_modeler.server import Server
from adze_modeler.simulation import sim
from model import BLDCMotor


@sim.register('basic')
def basic_run(model, modelparams, simprams):
    m = model(**modelparams)
    return m()

@sim.register('cogging')
def cogging_torque(model, modelparams, simprams):
    return {"msg": "Cogging torque simulation works"}


if __name__=='__main__':

    # set the model for the simulation
    sim.set_model(BLDCMotor)
    model = Server(sim)
    model.run()
