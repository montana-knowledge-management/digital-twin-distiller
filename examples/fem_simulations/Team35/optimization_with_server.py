import json
from math import inf
from operator import itemgetter

import requests
from artap.algorithm_genetic import NSGAII
from artap.problem import Problem


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
        req = {"simulation": {"type": "default", "x": x}}
        url = "http://192.168.0.117:5000/process"
        result = requests.post(url=url, data=json.dumps(req))
        if result.status_code == 200:
            # result.json() -> {'res': {'f1': 0.00015563575447000003}}
            f1 = result.json().get("res").get("f1")
            print(f1)
            return [f1]
        else:
            print("Failed")
            return [inf]


if __name__ == "__main__":

    # Perform the optimization iterating over 100 times on 100 individuals.
    problem = CoilOptimizationProblem()
    algorithm = NSGAII(problem)
    algorithm.options["max_population_number"] = 20
    algorithm.options["max_population_size"] = 50
    try:
        algorithm.run()
        res = problem.individuals[-1]
        print(res.vector)
        print(res.costs)
    except KeyboardInterrupt:
        pass
