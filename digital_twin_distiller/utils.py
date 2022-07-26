import csv
import functools
import warnings
from itertools import tee, zip_longest
from math import atan2, fmod, pi, sqrt
from pathlib import Path
from statistics import fmean
from uuid import uuid4

import matplotlib.pyplot as plt
from numpy import linspace
from numpy.polynomial import Polynomial as P

__all__ = [
    "getID",
    "mirror_point",
    "mm2px",
    "mm2inch",
    "get_width_height",
    "setup_matplotlib",
    "inch2mm",
    "rms",
    "csv_write",
    "get_polyfit",
    "pairwise",
    "get_short_id",
    "purge_dir",
    "deprecated",
    "get_phi",
]


def getID():
    return int(uuid4())


def mirror_point(p1, p2, p3):
    """
    Mirror the p3 point on the p1 - p2 line.

    https://math.stackexchange.com/questions/2192124/how-to-find-an-equation-a-mirroring-point-on-2d-space-mark-by-a-line
    """
    p12 = p2 - p1
    p13 = p3 - p1
    H = p1 + ((p13 @ p12) / abs(p12) ** 2) * p12
    return H + (H - p3)


def mm2px(x):
    """
    Convert millimeters to pixels
    """
    return int(3.7795275591 * x)


def mm2inch(x):
    """
    Convert millimeters to inches
    """
    return 0.03937007874 * x


def get_width_height(type_="onehalf", aspect=(16, 10), unit="px"):
    """
    This function returns the width and the height of a figure in pixels based on the Elsevier
    reccomendations on figure sizes.
    https://www.elsevier.com/authors/policies-and-guidelines/artwork-and-media-instructions/artwork-sizing

    Parameters:
        type_: The type of the figure, can be "minimal",
               "single", "onehalf", "full", "double"
        aspect: This iterable specifies the aspect ratio of the figure.
    """

    types = {
        "minimal": 30,
        "single": 90,
        "onehalf": 140,
        "full": 190,
        "double": 190,
    }
    units = {"px": mm2px, "inch": mm2inch}

    if type_ not in types.keys():
        raise ValueError(
            f"Invalid keyword argument. Got {type_=}. "
            'Accepted values are: "minimal", "single",'
            '"onehalf", "full", "double".'
        )

    scaley = aspect[0] / aspect[1]
    width = types[type_]
    height = width / scaley

    if unit == "mm":
        return width, height
    else:
        return units[unit](width), units[unit](height)


# TODO: ki kéne tenni sablonba
def setup_matplotlib():
    """
    Sets up a matplotlib template.
    """
    plt.style.use(["default", "seaborn-bright"])
    w, h = get_width_height(type_="onehalf", aspect=(16, 9), unit="inch")
    plt.rcParams["figure.figsize"] = w, h
    plt.rcParams["lines.linewidth"] = 1

    plt.rc("font", size=8)  # controls default text sizes
    plt.rc("axes", titlesize=8)  # fontsize of the axes title
    plt.rc("axes", labelsize=8)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=8)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=8)  # fontsize of the tick labels
    plt.rc("legend", fontsize=8)  # legend fontsize
    plt.rc("figure", titlesize=8)  # fontsize of the figure title

    # plt.grid(b=True, which="major", color="#666666", linestyle="-", linewidth=0.8)
    # plt.grid(b=True, which="minor", color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
    # plt.minorticks_on()
    # plt.xlabel("")
    # plt.ylabel("")
    # plt.legend()
    # plt.savefig(dir_media / ".pdf", bbox_inches="tight")
    # plt.show()


def inch2mm(x):
    """
    Convert inches to millimeters.
    """
    return 25.4 * x


def rms(arr):
    """
    Get the root mean square value from a Sequence.
    """
    return sqrt(fmean(map(lambda xi: xi**2, arr)))


def csv_write(file, names, *args):
    """
    Write arrays into a csv file.

    Parameters:
        file: The name of the csv file.
        names: Sequence of strings that specifies the column names.
        args: 1D sequences

    Example:
        x = [0,1,2,3]
        y=[-1,2,-33,0]
        csv_write('data.csv', ['iteration', 'measurement'], x, y)
    """
    # TODO: stem check for file object
    file = Path(file)
    assert file.parent.exists(), f"There is no directory: {file.parent}"
    assert len(names) == len(args), (
        f"The number of names({len(names)}) and " f"the number of columns({len(args)}) are not equal."
    )
    with open(file, "w", newline="", encoding="UTF8") as f:
        w = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
        w.writerow(names)
        w.writerows(zip(*args))


def csv_read(file, dict_return=False):
    """
    Read data from csv files. The function doesn't check if the first row has the column names.

    Parameters:
        file: The name of the csv file.
        dict_return: If True, the function will return a dictionary with the column names, else
                     it will return the data.

    """
    file = Path(file)
    assert file.exists(), f"File does not exists: {file}"
    with open(file, encoding="UTF8") as f:
        r = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
        names = next(r)
        data = zip(*(li for li in r))
        if dict_return:
            return {ni: di for ni, di in zip(names, data)}
        else:
            return data


def get_polyfit(x, y, N=None, verbose=False):
    """Fits a polynomial to the measurement points, gives back the refined points."""
    assert len(x) == len(y)
    N = 2 * len(x) + 1 if N is None else int(N)
    x_fine = linspace(min(x), max(x), N)
    maxy_ref = max(y)
    rmsy_ref = rms(y)
    cases = []
    for order in range(2, 19):
        p = P.fit(x, y, order)
        y_fine = p(x_fine)
        maxy_fine = max(y_fine)
        rmsy_fine = rms(y_fine)
        d1 = (maxy_fine - maxy_ref) / maxy_ref * 100
        d2 = (rmsy_fine - rmsy_ref) / rmsy_ref * 100
        cases.append((d1, d2, order))

    cases.sort(key=lambda ci: abs(ci[0]))
    best_order = cases[0][-1]
    p = P.fit(x, y, best_order)
    y_best = p(x_fine)
    if verbose:
        print(f"Best: MAX: {cases[0][0]:.3f} % RMS: {cases[0][1]:.3f} % order: {best_order}")
    return x_fine, y_best


def pairwise(iterable, cycle=False):
    """
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    https://docs.python.org/3/library/itertools.html#itertools.pairwise

    If incudelast:
        pairwise('ABCDEFG') --> AB BC CD DE EF FG GA
    """

    a, b = tee(iterable)
    c = next(b, None)

    if cycle:
        return zip_longest(a, b, fillvalue=c)
    else:
        return zip(a, b)


def get_short_id(point, n: int = 6):
    """
    This function gives back the Node instances' ID in a more readable string format.
    This function is used during the various transformation functions that modify the geometry.

    :param point: A Node instance
    :type: Node
    :param n: The number of digits that is going to be returned from the Nodes long ID.
    :type: int
    """
    return hex(point.id)[-6:]


def purge_dir(location, force=False):
    """
    Delete ALL FILES AND DIRECTORIES under location including location itself.
    USE WITH CAUTION !
    """

    # get the location path
    # if its pointing to a file then the parent directory will be the basis.
    location = Path(location)
    if location.suffix:
        if force:
            location = location.parent
        else:
            raise RuntimeError(
                "Location is pointing to a file. You can delete it alongside with it's parent directory"
                "by using the force=True flag. BE CAREFUL! It will delete everything that is located"
                " under the files parent directory!"
            )
    if not location.exists():
        return False

    # print('Deleting everything under:', location)
    dirs = []
    # iterating over everything that rglob finds
    for item in location.rglob("**/*"):
        # if item is a file we delete it
        if item.suffix:
            item.unlink()
        else:
            # if item is a directory we append it to a list
            dirs.append(item)

    # we sort the directory list based on their length
    # Note: Not the string length matters but the part count
    dirs.sort(key=lambda di: len(di.parts), reverse=True)
    # we delete the directories from the most nested one to the top level
    for dir in dirs:
        dir.rmdir()

    # lastly we delete location itself
    location.rmdir()


def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.
    Credit: https://stackoverflow.com/questions/2536307/decorators-in-the-python-standard-lib-deprecated-specifically
    """

    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.simplefilter("always", DeprecationWarning)  # turn off filter
        warnings.warn(f"Call to deprecated function {func.__name__}.", category=DeprecationWarning, stacklevel=2)
        warnings.simplefilter("default", DeprecationWarning)  # reset filter
        return func(*args, **kwargs)

    return new_func


def get_phi(x, y, tol=1e-5):
    """
    This function calculates the angle between a point and the x axis.
    :return [0, 360]
    """
    phi = atan2(y, x) * 180.0 / pi

    if phi < -tol:
        phi += 360.0

    phi = fmod(phi + tol, 360.0) - tol

    return phi
