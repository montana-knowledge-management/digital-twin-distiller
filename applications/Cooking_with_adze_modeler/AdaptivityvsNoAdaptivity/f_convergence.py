import math
import operator
from adze_modeler.metadata import Agros2DMetadata
from adze_modeler.metadata import FemmMetadata
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from pathlib import Path
from random import choice
from random import seed
from random import uniform

from problem_solver import build


def evaluate(X, a_tol):
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
    agros_metadata.nb_refinements = 0  # should be 0 otherwise adaptivity won't work !
    agros_metadata.polyorder = 2
    agros_metadata.adaptivity = "hp-adaptivity"
    agros_metadata.adaptivity_steps = 100
    agros_metadata.adaptivity_tol = a_tol
    platform_agros = Agros2D(agros_metadata)

    # platform = platform_femm
    platform = platform_agros

    if any(map(lambda i: i < 5, X[6:12])):
        print(X)
        raise ValueError()
    F1 = math.inf
    F2 = math.inf
    F3 = math.inf

    if not (res := build(platform, X, a_tol, cleanup=False)):
        return [math.inf, math.inf, math.inf]

    Bz = [pointvalue[2] for pointvalue in res["Bz"]]  # [x, y, Bz(x, y)]
    Br = [pointvalue[2] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
    xi = [pointvalue[0] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
    yi = [pointvalue[1] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
    nb_nodes = int(res["nodes"])

    # Calculate F1
    B0 = 2e-3
    F1 = max(map(lambda Bz_i: abs(Bz_i - B0), Bz))

    # Calculate F2
    perturbation = 0.5  # mm
    Xp = [xi + perturbation for xi in X]
    if not (resp := build(platform, Xp, a_tol)):
        return [math.inf, math.inf, math.inf]

    Bzp = [pointvalue[2] for pointvalue in resp["Bz"]]
    Brp = [pointvalue[2] for pointvalue in resp["Br"]]

    Xn = [xi - perturbation for xi in X]
    if not (resn := build(platform, Xn, a_tol)):
        return [math.inf, math.inf, math.inf]

    Bzn = [pointvalue[2] for pointvalue in resn["Bz"]]
    Brn = [pointvalue[2] for pointvalue in resn["Br"]]

    deltaBpz = map(operator.abs, map(operator.sub, Bzp, Bz))
    deltaBpr = map(operator.abs, map(operator.sub, Brp, Br))
    deltaBp = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBpz, deltaBpr))

    deltaBnz = map(operator.abs, map(operator.sub, Bzn, Bz))
    deltaBnr = map(operator.abs, map(operator.sub, Brn, Br))
    deltaBn = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBnz, deltaBnr))

    F2 = max(map(operator.add, deltaBp, deltaBn))

    # Calcukate F3
    F3 = sum(X)

    with open(Path(__file__).parent / "agros_adaptivity.csv", "a+") as f:
        """
        platform, F1, F2, F3, nodes, r0, r1, r2, r3, ..., r19
        """
        record = [res["platform"]]
        record.extend([F1, F2, F3])
        record.append(nb_nodes)
        if platform.metadata.compatible_platform == "agros2d":
            record.append(agros_metadata.adaptivity_tol)
        else:
            record.append(a_tol)

        f.write(",".join([str(i) for i in record]))
        f.write("\n")

    return [F1, F2, F3, nb_nodes]


X = [
    32.331913124436305,
    2.22552700591068,
    14.476436600086844,
    11.937326169292314,
    37.087089494036604,
    34.158274883722655,
    45.14808054671804,
    8.912247468323727,
    23.98648188583717,
    6.340874874713165,
    14.83870886616215,
    27.740987964651307,
    6.194118635773863,
    13.947694280899182,
    34.244799700078545,
    29.522366627144752,
    14.919827991831351,
    31.51695577441589,
    40.662092377213504,
    1.3184392242249898,
    40.485143339807586,
    35.208830354423114,
    17.6722753093816,
    8.618495490777295,
]


for i in range(10):
    tol = uniform(0, 1e-3)
    print(i + 1, tol, evaluate(X, tol))
