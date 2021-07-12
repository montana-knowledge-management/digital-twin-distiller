from random import uniform

from artap.problem import Problem
from artap.algorithm_genetic import NSGAII
from pathlib import Path
from adze_modeler.modelpiece import ModelPiece
from numpy import linspace
from problem_solver import generate_snapshot
from math import inf
basepath = Path(__file__).parent

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
        if s.execute(timeout=10) is not None:
            return s.retrive_results()['Fx']
        else:
            raise ValueError("Failed computation")


    def evaluate(self, individual):
        # X = individual.vector
        X=individual.copy()

        snapshots = [generate_snapshot(self.stator, self.rotor, X, di) for di in linspace(0, 50, 3)]
        try:
            F = [self.execute_snapshot(s) for s in snapshots]
        except ValueError:
            print("Gotcha.", F)
            return inf

        else:
            return F



if __name__=='__main__':
    bounds = ((2, 12),
              (0.1, 1),
              (6, 20.44),
              (0.1, 1),
              (24.44, 35.5),
              (0.1, 1),
              (30, 42.88),
              (0.1, 1),
              (46.88, 60),
              (0.1, 1),
              (52, 65.32),
              (0.1, 1))

    X = [uniform(lower*0.01, upper) for lower, upper in bounds]

    p = CoggingTorqueOptimizationProblem()
    p.set()
    p.evaluate(X)