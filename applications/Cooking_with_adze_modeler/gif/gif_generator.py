import operator
from adze_modeler.geometry import Geometry
from adze_modeler.modelpiece import ModelPiece
from pathlib import Path

from numpy import array
from numpy import unique


path_base = Path(__file__).parent.parent
current_dir = Path(__file__).parent
path_export_base = Path(__file__).parent / "media"


def build(X, nb):

    mp_bound = ModelPiece("bound")
    mp_bound.load_piece_from_svg(path_base / "ArtapOptimization" / "problem_boundary.svg")

    mp_control = ModelPiece("core")
    mp_control.load_piece_from_svg(path_base / "ArtapOptimization" / "core.svg")
    mp_control.put(0, -5)

    mp_coil = ModelPiece("coil")
    mp_coil.load_piece_from_svg(path_base / "ArtapOptimization" / "coil.svg")

    N = len(X)
    h = 1.5
    geom = Geometry()
    geom.merge_geometry(mp_bound.geom)
    geom.merge_geometry(mp_control.geom)

    for i, ri in enumerate(X):
        coil = mp_coil.spawn()
        offsetz = i * h - N * h / 2
        coil.put(ri, offsetz)
        geom.merge_geometry(coil.geom)

    geom.generate_intersections()
    geom.export_svg(path_export_base / f"nb-{nb:02d}.svg")


def get_line_from_file(filename):
    with open(filename) as f:
        yield from f


def get_processed_line_asym(filename):
    for line_i in get_line_from_file(filename):
        platform, *line_i = line_i.strip().split(",")
        f1, f2, f3, nodes, *r = (float(ri) for ri in line_i)
        yield (f1, f2, f3, *r)


def get_processed_line_sym(filename):
    for line_i in get_line_from_file(filename):
        platform, *line_i = line_i.strip().split(",")
        f1, f2, f3, nodes, *r = (float(ri) for ri in line_i)
        r = list(reversed(r))
        r.extend(reversed(r))
        yield (f1, f2, f3, *r)


sym_data_generator = get_processed_line_sym(Path(__file__).parent / "statistics_sym.csv")
asym_data_generator = get_processed_line_asym(Path(__file__).parent / "statistics_asym.csv")

data_sym = array([record for record in sym_data_generator])
data_asym = array([record for record in asym_data_generator])
data_asym = unique(data_asym, axis=0)
data_sym = unique(data_sym, axis=0)
idxF1 = 0
idxF2 = 1
idxF3 = 2

data_sym = array(sorted(data_sym, key=operator.itemgetter(idxF3)))
data_asym = array(sorted(data_asym, key=operator.itemgetter(idxF3)))

for i, xi in enumerate(data_asym[::100]):
    build(xi[3:], i)
    # convert '*.*' converted_%04d.png
