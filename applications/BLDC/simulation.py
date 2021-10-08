from adze_modeler.server import Server
from adze_modeler.simulation import sim
from model import BLDCMotor
from multiprocessing import Pool
from time import perf_counter
from numpy import linspace

def execute_model(model: BLDCMotor):
    t0 = perf_counter()
    result = model(timeout=2000, cleanup=False)
    t1 = perf_counter()
    result["Torque"] *= 8
    return result


@sim.register('basic')
def basic_run(model, modelparams, simprams, miscparams):
    m = model(**modelparams)
    result = m(cleanup=miscparams['cleanup'])
    result['Torque'] *= 8
    return result

@sim.register('cogging')
def cogging_torque(model, modelparams, simprams, miscparams):
    theta0 = simprams.get('theta0', 0)
    theta1 = simprams.get('theta1', 360/24)
    nsteps = simprams.get('nsteps', 12)
    theta = linspace(theta0, theta1, nsteps)
    models = []
    for theta_i in theta:
        modelparams['rotorangle'] = theta_i
        models.append(model(**modelparams))

    with Pool(processes=miscparams['processes']) as pool:
        results = pool.map(execute_model, models)

    return results


if __name__=='__main__':

    # set the model for the simulation
    sim.set_model(BLDCMotor)
    model = Server(sim)
    model.run()
