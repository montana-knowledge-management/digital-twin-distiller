from adze_modeler.modelpiece import ModelPiece
from pathlib import Path
from random import uniform
from uuid import uuid4

import matplotlib.pyplot as plt
from numpy import linspace
from problem_solver import generate_snapshot

basepath = Path(__file__).parent
exportpath = basepath / "snapshots"


def execute_snapshot(self, s):
    s.export()
    if s.execute(timeout=20) is not None:
        return s.retrive_results()["Fx"]
    else:
        raise ValueError("Failed computation")


def evaluate(X):
    computation_dir = exportpath / str(uuid4())
    computation_dir.mkdir(exist_ok=True)

    stator = ModelPiece("stator")
    stator.load_piece_from_dxf(basepath / "stator_stripped.dxf")
    stator.put(0, 0)
    rotor = ModelPiece("rotor")
    rotor.load_piece_from_dxf(basepath / "rotor.dxf")

    displacement = linspace(0, 50, 61)
    snapshots = [generate_snapshot(stator, rotor, X, di, export_loc=computation_dir) for di in displacement]

    T = []
    for si in snapshots:
        si.export()
        si.execute()
        T.append(si.retrive_results()["Fx"])

    T = [Ti * 300 / 1000 * 14 for Ti in T]
    plt.figure()
    plt.plot(displacement, T, "b-o")
    plt.xlabel("Displacement [mm]")
    plt.ylabel("Torque [Nm]")
    plt.grid()
    plt.savefig(basepath / "media" / "torque.png", dpi=550, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":

    bounds = (
        (2, 12),
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
        (0.1, 1),
    )

    X = [uniform(lower, upper) for lower, upper in bounds]
    for i in range(1, 13, 2):
        X[i] = 0.75

    evaluate(X)
