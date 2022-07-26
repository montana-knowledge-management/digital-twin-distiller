import json
from multiprocessing import Pool

from model import PriusMotor
from numpy import linspace

from digital_twin_distiller.encapsulator import Encapsulator
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import sim


def execute_model(model: PriusMotor):
    return model(timeout=2000, cleanup=True).get("Torque", 0.0) * 8


@sim.register("default")
def default_simulation(model, modelparams, simparams, miscparams):
    return "Hello World!"


@sim.register("cogging")
def cogging_calculation(model, modelparams, simparams, miscparams):
    alpha0 = simparams["alpha0"]
    alpha1 = simparams["alpha1"]
    nsteps = simparams["nsteps"]

    alpha = linspace(alpha0, alpha1, nsteps)
    models = [model(rotorangle=ai) for ai in alpha]

    with Pool() as pool:
        res = pool.map(execute_model, models)

    result = {"alpha": list(alpha), "T": list(res)}

    with open(ModelDir.DATA / f"cogging_torque.json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=True)

    return result


@sim.register("locked_rotor")
def locked_rotor(model, modelparams, simparams, miscparams):
    alpha0 = simparams["alpha0"]
    alpha1 = simparams["alpha1"]
    nsteps = simparams["nsteps"]
    I0 = simparams["I0"]

    alpha = linspace(alpha0, alpha1, nsteps)
    models = [model(I0=I0, alpha=ai) for ai in alpha]
    with Pool() as pool:
        T = pool.map(execute_model, models)

    result = {"I0": I0, "alpha": list(alpha), "T": list(T)}

    with open(ModelDir.DATA / f"locked_rotor.json", "w", encoding="utf-8") as f:
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
