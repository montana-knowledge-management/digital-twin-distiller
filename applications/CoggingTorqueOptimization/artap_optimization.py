from random import uniform
import operator
from artap.problem import Problem
from artap.algorithm_genetic import NSGAII
from pathlib import Path
from adze_modeler.modelpiece import ModelPiece
from numpy import linspace
from problem_solver import generate_snapshot
from math import inf
import matplotlib.pyplot as plt
import multiprocessing
from time import perf_counter
from shutil import rmtree
from uuid import uuid4

basepath = Path(__file__).parent
exportpath = basepath / "snapshots"


class CoggingTorqueOptimizationProblem(Problem):

    def set(self):
        self.name = 'Biobjective Test Problem'

        self.parameters = [{'name': "c11x", "bounds": [2, 12]},
                           {'name': "c11y", "bounds": [0.1, 1]},
                           {'name': "c21x", "bounds": [6, 20.44]},
                           {'name': "c21y", "bounds": [0.1, 1]},
                           {'name': "c12x", "bounds": [24.44, 35.5]},
                           {'name': "c12y", "bounds": [0.1, 1]},
                           {'name': "c22x", "bounds": [30, 42.88]},
                           {'name': "c22y", "bounds": [0.1, 1]},
                           {'name': "c13x", "bounds": [46.88, 60]},
                           {'name': "c13y", "bounds": [0.1, 1]},
                           {'name': "c23x", "bounds": [52, 65.32]},
                           {'name': "c23y", "bounds": [0.1, 1]}

                           ]

        self.costs = [{'name': 'f_1', 'criteria': 'minimize'}]

        self.stator = ModelPiece("stator")
        self.stator.load_piece_from_dxf(basepath / "stator_stripped.dxf")
        self.stator.put(0, 0)
        self.rotor = ModelPiece("rotor")
        self.rotor.load_piece_from_dxf(basepath / "rotor.dxf")

    def execute_snapshot(self, s):
        s.export()
        if s.execute(timeout=20) is not None:
            return s.retrive_results()['Fx']
        else:
            raise ValueError("Failed computation")


    def evaluate(self, individual, export=False):
        computation_dir = exportpath / str(uuid4())
        computation_dir.mkdir(exist_ok=True)

        X = individual.vector
        # X=individual.copy()
        displacement = linspace(0, 50, 25)
        snapshots = [generate_snapshot(self.stator, self.rotor, X, di, export_loc=computation_dir) for di in displacement]
        try:
            # t0 = perf_counter()
            with multiprocessing.Pool(processes=2) as pool:
                Fx = pool.map(self.execute_snapshot, snapshots)

            # t1 = perf_counter()
            # print(t1-t0)
            T = [Fi*300/1000*14 for Fi in Fx]
            F1 = sum(map(lambda Ti: Ti * Ti, T))
            with open(Path(__file__).parent / "statistics.csv", "a+") as f:
                record = [F1]
                record.extend(X)
                f.write(','.join([str(ri) for ri in record]))
                f.write('\n')

            if export:
                with open(Path(__file__).parent / "torque.csv", "w") as f:
                    for di, Ti in zip(displacement, T):
                        f.write(f'{di}, {Ti}\n')



        except ValueError:
            print("Gotcha.")
            return [inf]

        finally:
            rmtree(computation_dir)

        return [F1]



if __name__=='__main__':
    # bounds = ((2, 12),
    #           (0.1, 1),
    #           (6, 20.44),
    #           (0.1, 1),
    #           (24.44, 35.5),
    #           (0.1, 1),
    #           (30, 42.88),
    #           (0.1, 1),
    #           (46.88, 60),
    #           (0.1, 1),
    #           (52, 65.32),
    #           (0.1, 1))
    #
    # X = [uniform(lower, upper) for lower, upper in bounds]
    # for i in range(1, 13, 2):
    #     X[i] = 0.75

    p = CoggingTorqueOptimizationProblem()

    algorithm = NSGAII(p)
    algorithm.options['max_population_number'] = 3
    algorithm.options['max_population_size'] = 3
    try:
        algorithm.run()
        res = p.individuals[-1]
        print(res.vector)
        print(res.costs)
    except KeyboardInterrupt:
        pass

    with open(Path(__file__).parent / "pareto_front.csv", "w") as f:
        for ind in p.individuals:
            record = ind.costs.copy()
            record.extend(ind.vector.copy())
            f.write(','.join([str(i) for i in record]))
            f.write('\n')


    p.evaluate(p.individuals[-1], export=True)

    dstart = []
    Tstart = []
    with open(Path(__file__).parent / "torque_start.csv", "r") as f:
        for line in f.readlines():
            di, Ti = line.strip().split(',')
            dstart.append(float(di))
            Tstart.append(float(Ti))

    dbest = []
    Tbest = []
    with open(Path(__file__).parent / "torque.csv", "r") as f:
        for line in f.readlines():
            di, Ti = line.strip().split(',')
            dbest.append(float(di))
            Tbest.append(float(Ti))


    plt.figure()
    plt.plot(dstart, Tstart, 'b-o', label="Start")
    plt.plot(dbest, Tbest, 'r-o', label="End")
    plt.legend()
    plt.grid()
    plt.savefig(Path(__file__).parent / "media" / "torque.png")
    plt.show()