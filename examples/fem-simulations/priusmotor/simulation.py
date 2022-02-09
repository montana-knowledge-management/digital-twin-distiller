import json
from multiprocessing import Pool
from itertools import product
from numpy import linspace

from digital_twin_distiller.encapsulator import Encapsulator
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import sim
from model import PriusMotor

def execute_model(model: PriusMotor):
    return model(timeout=2000, cleanup=True).get("Torque", 0.0) * -8

@sim.register('default')
def default_simulation(model, modelparams, simparams, miscparams):
    return "Hello World!"

@sim.register('cogging')
def cogging_calculation(model, modelparams, simparams, miscparams):
    range_a0 = simparams['range_a0']
    range_a1 = simparams['range_a1']
    nsteps_a = simparams['nsteps_a']

    range_b0 = simparams['range_b0']
    range_b1 = simparams['range_b1']
    nsteps_b = simparams['nsteps_b']

    range_c0 = simparams['range_c0']
    range_c1 = simparams['range_c1']
    nsteps_c = simparams['nsteps_c']

    range_a = linspace(range_a0, range_a1, nsteps_a)
    range_b = linspace(range_b0, range_b1, nsteps_b)
    range_c = linspace(range_c0, range_c1, nsteps_c)

    prod = list(product(range_a, range_b, range_c))

    models = [model(earheight=ai, aslheight=bi, rotorangle=ci) for ai in range_a for bi in range_b for ci in range_c]
    with Pool() as pool:
        res = pool.map(execute_model, models)

    result = {'Torque': list(res)}

    with open(ModelDir.DATA / f'cogging_torque.json', 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2, ensure_ascii=True)

    return result

@sim.register('locked_rotor')
def locked_rotor(model, modelparams, simparams, miscparams):
    range_a0 = simparams['range_a0']
    range_a1 = simparams['range_a1']
    nsteps_a = simparams['nsteps_a']

    range_b0 = simparams['range_b0']
    range_b1 = simparams['range_b1']
    nsteps_b = simparams['nsteps_b']

    range_c0 = simparams['range_c0']
    range_c1 = simparams['range_c1']
    nsteps_c = simparams['nsteps_c']

    range_a = linspace(range_a0, range_a1, nsteps_a)
    range_b = linspace(range_b0, range_b1, nsteps_b)
    range_c = linspace(range_c0, range_c1, nsteps_c)

    I0 = simparams["I0"]

    prod = list(product(range_a, range_b, range_c))

    models = [model(I0=I0, earheight=ai, aslheight=bi, rotorangle=ci) for ai in range_a for bi in range_b for ci in range_c]
    with Pool() as pool:
        res = pool.map(execute_model, models)

    result = {'Torque': list(res)}

    with open(ModelDir.DATA / f'locked_rotor.json', 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2, ensure_ascii=True)

    return result

if __name__ == "__main__":
    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(PriusMotor)

    model = Encapsulator(sim)
    model.port = 8080
    model.build_docs()
    model.run()
