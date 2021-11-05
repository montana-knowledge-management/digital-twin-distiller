import math
import multiprocessing
import operator
from functools import partial
from pathlib import Path
from random import choice, uniform

from artap.algorithm_genetic import NSGAII
from artap.problem import Problem
from symmetric_model import SymmetircModel

from adze_modeler.metadata import Agros2DMetadata, FemmMetadata
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm

DIR_BASE = Path(__file__).parent
DIR_DATA = DIR_BASE / "data"


def execute_model(model_i):
    return model_i()


class CoilOptimizationProblem(Problem):
    def set(self):
        self.name = "Biobjective Test Problem"

        self.parameters = [
            {"name": "r0", "bounds": [5.5, 20]},
            {"name": "r1", "bounds": [5.5, 20]},
            {"name": "r2", "bounds": [5.5, 20]},
            {"name": "r3", "bounds": [5.5, 20]},
            {"name": "r4", "bounds": [5.5, 20]},
            {"name": "r5", "bounds": [5.5, 20]},
            {"name": "r6", "bounds": [5.5, 20]},
            {"name": "r7", "bounds": [5.5, 20]},
            {"name": "r8", "bounds": [5.5, 20]},
            {"name": "r9", "bounds": [5.5, 20]},
        ]

        self.costs = [
            {"name": "f_1", "criteria": "minimize"},
            {"name": "f_2", "criteria": "minimize"},
            {"name": "f_3", "criteria": "minimize"},
        ]

    def evaluate(self, individual):
        perturbation = 0.5  # mm
        X = individual.vector
        Xp = [xi + perturbation for xi in X]
        Xn = [xi - perturbation for xi in X]
        models = [SymmetircModel(X), SymmetircModel(Xn), SymmetircModel(Xp)]

        F1 = math.inf
        F2 = math.inf
        F3 = math.inf
        try:

            with multiprocessing.Pool(processes=3) as pool:
                res, resn, resp = pool.map(execute_model, models)

            Bz = [pointvalue[2] for pointvalue in res["Bz"]]  # [x, y, Bz(x, y)]
            Br = [pointvalue[2] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
            xi = [pointvalue[0] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
            yi = [pointvalue[1] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
            nb_nodes = res["nodes"]

            # Calculate F1
            B0 = 2e-3
            F1 = max(map(lambda Bz_i: abs(Bz_i - B0), Bz))

            # Calculate F2

            Bzp = [pointvalue[2] for pointvalue in resp["Bz"]]
            Brp = [pointvalue[2] for pointvalue in resp["Br"]]

            Bzn = [pointvalue[2] for pointvalue in resn["Bz"]]
            Brn = [pointvalue[2] for pointvalue in resn["Br"]]

            deltaBpz = map(operator.abs, map(operator.sub, Bzp, Bz))
            deltaBpr = map(operator.abs, map(operator.sub, Brp, Br))
            deltaBp = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBpz, deltaBpr))

            deltaBnz = map(operator.abs, map(operator.sub, Bzn, Bz))
            deltaBnr = map(operator.abs, map(operator.sub, Brn, Br))
            deltaBn = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBnz, deltaBnr))

            F2 = max(map(operator.add, deltaBp, deltaBn))

            # Calculate F3
            F3 = sum(X)

            print(F1, F2, F3)
            return [F1, F2, F3]
        except Exception as e:
            return [math.inf, math.inf, math.inf]


if __name__ == "__main__":

    # Perform the optimization iterating over 100 times on 100 individuals.
    problem = CoilOptimizationProblem()
    algorithm = NSGAII(problem)
    algorithm.options["max_population_number"] = 100
    algorithm.options["max_population_size"] = 100
    try:
        algorithm.run()
        res = problem.individuals[-1]
        print(res.vector)
        print(res.costs)
    except KeyboardInterrupt:
        pass

    with open(DIR_DATA / "pareto_front_symmetric.csv", "w") as f:
        for ind in problem.individuals:
            record = ind.costs.copy()
            record.extend(ind.vector.copy())
            f.write(",".join([str(i) for i in record]))
            f.write("\n")
