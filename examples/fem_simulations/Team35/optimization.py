from math import inf
from operator import itemgetter

from artap.algorithm_genetic import NSGAII
from artap.problem import Problem
from model import DistributedWinding


class CoilOptimizationProblem(Problem):
    def set(self):
        self.name = "Team35 Test Problem"

        self.parameters = [
            {"name": "r0", "bounds": [5, 20]},
            {"name": "r1", "bounds": [5, 20]},
            {"name": "r2", "bounds": [5, 20]},
            {"name": "r3", "bounds": [5, 20]},
            {"name": "r4", "bounds": [5, 20]},
            {"name": "r5", "bounds": [5, 20]},
            {"name": "r6", "bounds": [5, 20]},
            {"name": "r7", "bounds": [5, 20]},
            {"name": "r8", "bounds": [5, 20]},
            {"name": "r9", "bounds": [5, 20]},
        ]

        self.costs = [{"name": "f_1", "criteria": "minimize"}]

    def evaluate(self, individual):
        x = individual.vector

        x1 = [round(xi, 2) for xi in x]
        print("called with", x1, end=" ")
        assert len(x) == 10

        try:
            model = DistributedWinding(x)
            res = model(devmode=False, timeout=30, cleanup=True)

            # Br = map(op.itemgetter(2), res["Br"])
            Bz = map(itemgetter(2), res["Bz"])

            B0 = 2e-3
            f1 = max(map(lambda Bz_i: abs(Bz_i - B0), Bz))

            print("DONE")
            return [f1]
        except:
            print("FAILED")
            return [inf]


if __name__ == "__main__":

    # Perform the optimization iterating over 100 times on 100 individuals.
    problem = CoilOptimizationProblem()
    algorithm = NSGAII(problem)
    algorithm.options["max_population_number"] = 10
    algorithm.options["max_population_size"] = 10
    try:
        algorithm.run()
        res = problem.individuals[-1]
        print(res.vector)
        print(res.costs)
    except KeyboardInterrupt:
        pass
