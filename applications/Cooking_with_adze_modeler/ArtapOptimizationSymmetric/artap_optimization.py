import operator
from pathlib import Path
from artap.results import Results
import math
from artap.problem import Problem
from artap.algorithm_genetic import NSGAII

from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from problem_solver import build
from adze_modeler.metadata import Agros2DMetadata
from adze_modeler.metadata import FemmMetadata
from random import choice, uniform

current_dir = Path(__file__).parent

class CoilOptimizationProblem(Problem):

    def set(self):
        self.name = 'Biobjective Test Problem'

        self.parameters = [{'name': 'r0', 'bounds': [5.5, 50]},
                           {'name': 'r1', 'bounds': [5.5, 50]},
                           {'name': 'r2', 'bounds': [5.5, 50]},
                           {'name': 'r3', 'bounds': [5.5, 50]},
                           {'name': 'r4', 'bounds': [1, 50]},
                           {'name': 'r5', 'bounds': [1, 50]},
                           {'name': 'r6', 'bounds': [1, 50]},
                           {'name': 'r7', 'bounds': [1, 50]},
                           {'name': 'r8', 'bounds': [1, 50]},
                           {'name': 'r9', 'bounds': [1, 50]}
                           ]

        self.costs = [{'name': 'f_1', 'criteria': 'minimize'},
                      {'name': 'f_2', 'criteria': 'minimize'},
                      {'name': 'f_3', 'criteria': 'minimize'}]

    def evaluate(self, individual):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "axisymmetric"
        femm_metadata.file_script_name = "femm_solver_script"
        femm_metadata.file_metrics_name = "femm_solution.csv"
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = False
        platform_femm = Femm(femm_metadata)

        agros_metadata = Agros2DMetadata()
        agros_metadata.file_script_name = "agros_solver_script"
        agros_metadata.file_metrics_name = "agros_solution.csv"
        agros_metadata.problem_type = "magnetic"
        agros_metadata.coordinate_type = "axisymmetric"
        agros_metadata.analysis_type = "steadystate"
        agros_metadata.unit = 1e-3
        agros_metadata.nb_refinements = 0
        agros_metadata.polyorder = 2
        agros_metadata.adaptivity = "hp-adaptivity"
        agros_metadata.adaptivity_tol = 1
        platform_agros = Agros2D(agros_metadata)

        platform = platform_agros

        X = individual.vector
        F1 = math.inf
        F2 = math.inf
        F3 = math.inf

        if not (res := build(platform, X)):
            return [math.inf, math.inf, math.inf]

        Bz = [pointvalue[2] for pointvalue in res['Bz']] # [x, y, Bz(x, y)]
        Br = [pointvalue[2] for pointvalue in res['Br']] # [x, y, Br(x, y)]
        xi = [pointvalue[0] for pointvalue in res['Br']] # [x, y, Br(x, y)]
        yi = [pointvalue[1] for pointvalue in res['Br']] # [x, y, Br(x, y)]
        nb_nodes = res['nodes']

        # Calculate F1
        B0 = 2e-3
        F1 = max(map(lambda Bz_i: abs(Bz_i - B0), Bz))

        # Calculate F2
        perturbation = 0.5 # mm
        Xp = [xi + perturbation for xi in X]
        if not (resp := build(platform, Xp)):
            return [math.inf, math.inf, math.inf]

        Bzp = [pointvalue[2] for pointvalue in resp['Bz']]
        Brp = [pointvalue[2] for pointvalue in resp['Br']]

        Xn = [xi - perturbation for xi in X]
        if not (resn := build(platform, Xn)):
            return [math.inf, math.inf, math.inf]

        Bzn = [pointvalue[2] for pointvalue in resn['Bz']]
        Brn = [pointvalue[2] for pointvalue in resn['Br']]

        deltaBpz = map(operator.abs, map(operator.sub, Bzp, Bz))
        deltaBpr = map(operator.abs, map(operator.sub, Brp, Br))
        deltaBp = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBpz, deltaBpr))


        deltaBnz = map(operator.abs, map(operator.sub, Bzn, Bz))
        deltaBnr = map(operator.abs, map(operator.sub, Brn, Br))
        deltaBn = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBnz, deltaBnr))

        F2 = max(map(operator.add, deltaBp, deltaBn))

        # Calculate F3
        F3 = sum(X)

        with open(current_dir / 'statistics.csv', 'a+') as f:
            """
            platform, F1, F2, F3, nodes, r0, r1, r2, r3, ..., r19
            """
            record = [res['platform']]
            record.extend([F1, F2, F3])
            record.append(nb_nodes)
            record.extend(X)
            f.write(','.join([str(i) for i in record]))
            f.write('\n')

        # with open(current_dir / "trainingdata.csv", "a+") as f:
        #     """
        #     x, y, r0, r1, r2, ..., r19, Br, Bz
        #     """
        #     for x_i, y_i, Br_i, Bz_i in zip(xi, yi, Br, Bz):
        #         record = [x_i, y_i]
        #         record.extend(X)
        #         record.append(Br_i)
        #         record.append(Bz_i)
        #         f.write(','.join([str(i) for i in record]))
        #         f.write('\n')

        # print(F1, F2, F2)
        return [F1, F2, F3]

if __name__=='__main__':

    # Perform the optimization iterating over 100 times on 100 individuals.
    problem = CoilOptimizationProblem()
    algorithm = NSGAII(problem)
    algorithm.options['max_population_number'] = 2
    algorithm.options['max_population_size'] = 5
    try:
        algorithm.run()
        res = problem.individuals[-1]
        print(res.vector)
        print(res.costs)
    except KeyboardInterrupt:
        pass

    # with open(Path(__file__).parent / "pareto_front.csv", "w") as f:
    #     for ind in problem.individuals:
    #         record = ind.costs.copy()
    #         record.extend(ind.vector.copy())
    #         f.write(','.join([str(i) for i in record]))
    #         f.write('\n')


    # surrogate modelles megoldás