import itertools
from typing import Sequence

import numpy as np
from scipy.linalg import hankel, toeplitz

"""
This code was originally published by the following individuals for use with
Scilab:
    Copyright (C) 2012 - 2013 - Michael Baudin
    Copyright (C) 2012 - Maria Christopoulou
    Copyright (C) 2010 - 2011 - INRIA - Michael Baudin
    Copyright (C) 2009 - Yann Collette
    Copyright (C) 2009 - CEA - Jean-Marc Martinez
    website: forge.scilab.org/index.php/p/scidoe/sourcetree/master/macros
Credit goes to these individuals.

https://github.com/tirthajyoti/Design-of-experiment-Python
https://pythonhosted.org/pyDOE/

GSD Copyright (C) 2018 - Rickard Sjoegren
https://github.com/clicumu/pyDOE2/blob/master/pyDOE2

"""

"""

https://github.com/tirthajyoti/Design-of-experiment-Python/blob/fc1d00b9525e7e583153727a8979b9427122a3e4/pyDOE_corrected.py#L248

https://github.com/tisimst/pyDOE/blob/master/pyDOE/doe_plackett_burman.py

"""

__all__ = ["doe_fullfact", "doe_bbdesign", "doe_pbdesign", "doe_ccf"]


def doe_fullfact(levels: Sequence[int]):
    """
    Generate a general full-factorial design

    Parameters
    ----------
    levels : Sequence of integers
        An array of integers that indicate the number of levels of each input
        design factor.

    Returns
    -------
    mat : 2d-array
        The design matrix with coded levels 0 to k-1 for a k-level factor
    """

    return np.array(list(itertools.product(*(range(ni) for ni in levels))))


def ff2n(n: int = 1):
    """
    Create a 2-Level full-factorial design

    Parameters
    ----------
    n : int
        The number of factors in the design.

    Returns
    -------
    mat : 2d-array
        The design matrix with coded levels -1 and 1
    """
    return 2 * doe_fullfact([2] * n) - 1


def doe_bbdesign(n: int = 3, center=None):
    assert n >= 3, "Number of variables must be at least 3"

    repeat_center = lambda n, repeat: np.zeros((repeat, n))
    H_fact = ff2n(2)

    index = 0
    nb_lines = int((0.5 * n * (n - 1)) * H_fact.shape[0])
    H = repeat_center(n, nb_lines)

    for i in range(n - 1):
        for j in range(i + 1, n):
            index = index + 1
            H[
                max([0, (index - 1) * H_fact.shape[0]]) : index * H_fact.shape[0],
                i,
            ] = H_fact[:, 0]
            H[
                max([0, (index - 1) * H_fact.shape[0]]) : index * H_fact.shape[0],
                j,
            ] = H_fact[:, 1]

    if center is None:
        if n <= 16:
            points = [0, 0, 0, 3, 3, 6, 6, 6, 8, 9, 10, 12, 12, 13, 14, 15, 16]
            center = points[n]
        else:
            center = n

    H = np.c_[H.T, repeat_center(n, center).T].T

    return [list(li.astype(int)) for li in H]


def doe_pbdesign(n):
    """
    Plackett-Burman design
    """
    assert n > 0, "Number of factors must be a positive integer"
    keep = int(n)
    n = 4 * (int(n / 4) + 1)  # calculate the correct number of rows (multiple of 4)
    f, e = np.frexp([n, n / 12.0, n / 20.0])
    k = [idx for idx, val in enumerate(np.logical_and(f == 0.5, e > 0)) if val]

    assert isinstance(n, int) and k != [], "Invalid inputs. n must be a multiple of 4."

    k = k[0]
    e = e[k] - 1

    if k == 0:  # N = 1*2**e
        H = np.ones((1, 1))
    elif k == 1:  # N = 12*2**e
        H = np.vstack(
            (
                np.ones((1, 12)),
                np.hstack(
                    (
                        np.ones((11, 1)),
                        toeplitz(
                            [-1, -1, 1, -1, -1, -1, 1, 1, 1, -1, 1],
                            [-1, 1, -1, 1, 1, 1, -1, -1, -1, 1, -1],
                        ),
                    )
                ),
            )
        )
    elif k == 2:  # N = 20*2**e
        H = np.vstack(
            (
                np.ones((1, 20)),
                np.hstack(
                    (
                        np.ones((19, 1)),
                        hankel(
                            [
                                -1,
                                -1,
                                1,
                                1,
                                -1,
                                -1,
                                -1,
                                -1,
                                1,
                                -1,
                                1,
                                -1,
                                1,
                                1,
                                1,
                                1,
                                -1,
                                -1,
                                1,
                            ],
                            [
                                1,
                                -1,
                                -1,
                                1,
                                1,
                                -1,
                                -1,
                                -1,
                                -1,
                                1,
                                -1,
                                1,
                                -1,
                                1,
                                1,
                                1,
                                1,
                                -1,
                                -1,
                            ],
                        ),
                    )
                ),
            )
        )

    # Kronecker product construction
    for i in range(e):
        H = np.vstack((np.hstack((H, H)), np.hstack((H, -H))))

    # Reduce the size of the matrix as needed
    H = H[:, 1 : (keep + 1)]
    H = np.flipud(H)

    return [list(li.astype(int)) for li in H]


def doe_ccf(n):
    """
    Central Composite Design. In this design, the star points are at the center
    of each face of the factorial space. This variety requires 3 levels of each
    factor. Augmenting an existing factorial or resolution V design with
    appropriate star points can also produce this design.
    """

    H2 = np.zeros((2 * n, n))
    for i in range(n):
        H2[2 * i : 2 * i + 2, i] = [-1, 1]

    H1 = ff2n(n)

    ## redundant rows when (4, 4)
    center = (1, 0)
    C1 = np.zeros((center[0], n))
    C2 = np.zeros((center[1], n))
    H1 = np.r_[H1, C1]
    H2 = np.r_[H2, C2]
    H = np.r_[H1, H2]

    return [list(li.astype(int)) for li in H]


if __name__ == "__main__":
    import oapackage
    from numpy import array

    run_size = 18
    number_of_factors = 5
    factor_levels = 3
    strength = 5

    arrayclass = oapackage.arraydata_t(factor_levels, run_size, strength, number_of_factors)
    arr = array(arrayclass.create_root().getarray())
    print(arr + 1)
    # print(*doe_bbdesign(5, center=1), sep='\n')
    # print(*fullfact([3]*4), sep='\n')
    # print(*doe_pbdesign(4), sep='\n')
    # print(*dore_ccf(4), sep='\n')
