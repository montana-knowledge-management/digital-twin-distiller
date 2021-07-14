from pathlib import Path
import math
from artap.problem import Problem
from artap.algorithm_genetic import NSGAII
from artap.algorithm_nlopt import NLopt, LN_BOBYQA

from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from problem_solver import build
from adze_modeler.metadata import Agros2DMetadata
from adze_modeler.metadata import FemmMetadata

current_dir = Path(__file__).parent

class CoilOptimizationProblem(Problem):

    def set(self):
        self.name = 'Biobjective Test Problem'

        self.parameters = [{'name': 'r0', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r1', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r2', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r3', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r4', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r5', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r6', 'bounds': [5.5, 20], 'initial_value': 10},
                           {'name': 'r7', 'bounds': [5.5, 20], 'initial_value': 10},
                           {'name': 'r8', 'bounds': [5.5, 20], 'initial_value': 10},
                           {'name': 'r9', 'bounds': [5.5, 20], 'initial_value': 10},
                           {'name': 'r10', 'bounds': [5.5, 20], 'initial_value': 10},
                           {'name': 'r11', 'bounds': [5.5, 20], 'initial_value': 10},
                           {'name': 'r12', 'bounds': [5.5, 20], 'initial_value': 10},
                           {'name': 'r13', 'bounds': [5.5, 20], 'initial_value': 10},
                           {'name': 'r14', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r15', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r16', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r17', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r18', 'bounds': [1, 20], 'initial_value': 10},
                           {'name': 'r19', 'bounds': [1, 20], 'initial_value': 10}
                           ]

        self.costs = [{'name': 'f_1', 'criteria': 'minimize'}]

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
        try:

            X = individual.vector
            F1 = math.inf

            res = build(platform, X)

            Bz = [pointvalue[2] for pointvalue in res['Bz']] # [x, y, Bz(x, y)]
            Br = [pointvalue[2] for pointvalue in res['Br']] # [x, y, Br(x, y)]
            xi = [pointvalue[0] for pointvalue in res['Br']] # [x, y, Br(x, y)]
            yi = [pointvalue[1] for pointvalue in res['Br']] # [x, y, Br(x, y)]
            nb_nodes = res['nodes']

            # Calculate F1
            B0 = 2e-3
            nb_cols = 10
            nb_rows = 10
            F1 =0.0
            part1 = 0.0
            part2 = 0.0

            for n in range(1, nb_rows - 1):
                for m in range(1, nb_cols - 1):
                    part1 += math.sqrt((Bz[n * nb_cols + m] - Bz[(n - 1) * nb_cols + (m - 1)]) ** 2)
                    part1 += math.sqrt((Bz[n * nb_cols + m] - Bz[(n - 1) * nb_cols + m]) ** 2)
                    part1 += math.sqrt((Bz[n * nb_cols + m] - Bz[(n - 1) * nb_cols + (m + 1)]) ** 2)

                    part1 += math.sqrt((Bz[n * nb_cols + m] - Bz[n * nb_cols + (m - 1)]) ** 2)
                    part1 += math.sqrt((Bz[n * nb_cols + m] - Bz[n * nb_cols + m]) ** 2)
                    part1 += math.sqrt((Bz[n * nb_cols + m] - Bz[n * nb_cols + (m + 1)]) ** 2)

                    part1 += math.sqrt((Bz[n * nb_cols + m] - Bz[(n + 1) * nb_cols + (m - 1)]) ** 2)
                    part1 += math.sqrt((Bz[n * nb_cols + m] - Bz[(n + 1) * nb_cols + m]) ** 2)
                    part1 += math.sqrt((Bz[n * nb_cols + m] - Bz[(n + 1) * nb_cols + (m + 1)]) ** 2)

                    part2 += abs(Bz[n * nb_cols + m] - B0)

            F1 = part1 + part2

            with open(current_dir / 'statistics_bobyqa.csv', 'a+') as f:
                """
                platform, F1, nodes, r0, r1, r2, r3, ..., r19
                """
                record = [res['platform']]
                record.extend([F1])
                record.append(nb_nodes)
                record.extend(X)
                f.write(','.join([str(i) for i in record]))
                f.write('\n')

            print(F1)
            return [F1]

        except Exception as e:
            return [math.inf]

if __name__=='__main__':

    # Perform the optimization iterating over 100 times on 100 individuals.
    problem = CoilOptimizationProblem()
    # algorithm = NSGAII(problem)
    # algorithm.options['max_population_number'] = 100
    # algorithm.options['max_population_size'] = 100

    algorithm = NLopt(problem)
    algorithm.options['algorithm'] = LN_BOBYQA
    algorithm.options['n_iterations'] = 1000

    try:
        algorithm.run()
        res = problem.individuals[-1]
        print(res.vector)
        print(res.costs)
    except KeyboardInterrupt:
        pass

    with open(Path(__file__).parent / "pareto_front_bobyqa.csv", "w") as f:
        for ind in problem.individuals:
            record = ind.costs.copy()
            record.extend(ind.vector.copy())
            f.write(','.join([str(i) for i in record]))
            f.write('\n')