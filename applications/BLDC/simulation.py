from adze_modeler.modelpaths import ModelDir
from adze_modeler.server import Server
from adze_modeler.simulation import sim
from model import BLDCMotor
from multiprocessing import Pool
from time import perf_counter
from numpy import linspace
from random import normalvariate

def execute_model(model: BLDCMotor):
    t0 = perf_counter()
    result = model(timeout=2000, cleanup=False)
    t1 = perf_counter()
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

@sim.register('tol1')
def tol1(model, modelparams, simprams, miscparams):
    r1 = modelparams['r1']
    r2 = modelparams['r2']
    return {'Torque': normalvariate(r1, r2), "dummy": [0, 0, r1*r2]}

@sim.register('tol2')
def tol2(model, modelparams, simprams, miscparams):
    r1 = modelparams['r1']
    r2 = modelparams['r2']
    r = []
    for i in range(simprams['nsteps']):
        r.append({'Torque': normalvariate(r1, r2)+i, "dummy": [0, 0, i*r1 * r2]})
    return r



if __name__=='__main__':

    ModelDir.set_base(__file__)
    # set the model for the simulation
    sim.set_model(BLDCMotor)

    model = Server(sim)
    model.run()
