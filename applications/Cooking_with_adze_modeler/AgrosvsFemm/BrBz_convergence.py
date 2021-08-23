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


def evaluate(X, tol):
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
    agros_metadata.polyorder = 1
    agros_metadata.adaptivity = "hp-adaptivity"
    agros_metadata.adaptivity_steps = 100
    agros_metadata.adaptivity_tol = tol
    platform_agros = Agros2D(agros_metadata)

    platform = platform_femm
    # platform = platform_agros

    res = build(platform, X, tol, cleanup=True)

    Bz = [pointvalue[2] for pointvalue in res["Bz"]]  # [x, y, Bz(x, y)]
    Br = [pointvalue[2] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
    xi = [pointvalue[0] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
    yi = [pointvalue[1] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
    nb_nodes = int(res["nodes"])
    Bz_ = Bz[0]
    Br_ = Br[0]

    with open(Path(__file__).parent / "results.csv", "a+") as f:
        """
        platform, Br, Bz, nodes
        """
        record = [res["platform"]]
        record.extend([Br_, Bz_, nb_nodes])

        f.write(",".join([str(i) for i in record]))
        f.write("\n")

    return [Br_, Bz_, nb_nodes]

if __name__ == "__main__":

    X = [
        1.07757234384367,
        12.0810986190214,
        5.12071874756443,
        20.1720140666903,
        16.271996038583,
        5.52631142548928,
        25.7481784448498,
        9.35590470354692,
        45.9877661249235,
        7.40118928105421,
        10.4434974707302,
        38.4649413574651,
        16.6263411147178,
        7.76793550511193,
        8.42543646671979,
        22.6415022596837,
        46.928319345423,
        5.69151040344659,
        4.91528527435156,
        1.37343445767236,
    ]

    # for i in range(20):
    #     tol = uniform(1, 3)
    #     print(i + 1, tol, evaluate(X, tol))

    tol=0.1
    print(1, tol, evaluate(X, tol))