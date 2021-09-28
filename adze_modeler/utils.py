from uuid import uuid4
import matplotlib.pyplot as plt
from math import sqrt
from statistics import fmean
from pathlib import Path
import csv

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


def get_width_height(type_='onehalf', aspect=(16, 10), get_mm=False):
    """
    This function returns the width and the height of a figure in pixels based on the Elsevier
    reccomendations on figure sizes.
    https://www.elsevier.com/authors/policies-and-guidelines/artwork-and-media-instructions/artwork-sizing

    Parameters:
        type_: The type of the figure, can be "minimal",
               "single", "onehalf", "full", "double"
        aspect: This iterable specifies the aspect ratio of the figure.
    """
    width, height = 0, 0
    types = {'minimal': 30, 'single': 90, 'onehalf': 140, 'full': 190, 'double': 190}

    if type_ not in types.keys():
        raise ValueError(f'Invalid keyword argument. Got {type_=}. '
                         'Accepted values are: "minimal", "single",'
                         '"onehalf", "full", "double".')

    scaley = aspect[0] / aspect[1]
    width = types[type_]
    height = width / scaley

    if get_mm:
        return width, height
    else:
        return mm2px(width), mm2px(height)

def setup_matplotlib():
    plt.style.use(["default", "seaborn-bright"])
    w, h = get_width_height(type_='onehalf', aspect=(16, 9), get_mm=True)
    plt.rcParams["figure.figsize"] = mm2inch(w), mm2inch(h)
    plt.rcParams["lines.linewidth"] = 1

    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 14

    plt.rc("font", size=MEDIUM_SIZE)  # controls default text sizes
    plt.rc("axes", titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc("xtick", labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc("ytick", labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
    plt.rc("legend", fontsize=MEDIUM_SIZE)  # legend fontsize
    plt.rc("figure", titlesize=MEDIUM_SIZE)  # fontsize of the figure title

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
    assert len(names)==len(args), f"The number of names({len(names)}) and " \
                                  f"the number of columns({len(args)}) are not equal."
    with open(file, 'w', encoding='UTF8') as f:
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
    with open(file, 'r', encoding='UTF8') as f:
        r = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
        names = next(r)
        data = (li for li in r)
        if dict_return:
            return {ni: di for ni, di in zip(names, zip(*data))}
        else:
            return zip(*data)