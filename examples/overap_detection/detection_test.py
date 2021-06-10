# TODO: remove matplotlib after testing
import numpy as np
import matplotlib.pyplot as plt
from adze_modeler.geometry import Geometry

def get_intersetion_points(line_1, line_2, tol=1e-3):
    """

    https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

    Conditions:
    1. If r × s = 0 and (q − p) × r = 0, then the two lines are collinear.
    2. If r × s = 0 and (q − p) × r ≠ 0, then the two lines are parallel and non-intersecting.
    3. If r × s ≠ 0 and 0 ≤ t ≤ 1 and 0 ≤ u ≤ 1, the two line segments meet at the point p + t r = q + u s.
    4. Otherwise, the two line segments are not parallel but do not intersect.
    """
    x1 = line_1.start_pt.x
    y1 = line_1.start_pt.y
    x2 = line_1.end_pt.x
    y2 = line_1.end_pt.y

    x3 = line_2.start_pt.x
    y3 = line_2.start_pt.y
    x4 = line_2.end_pt.x
    y4 = line_2.end_pt.y

    p1 = None
    p2 = None

    p = np.array([x1, y1])
    r = np.array([x2 - x1, y2 - y1])

    q = np.array([x3, y3])
    s = np.array([x4 - x3, y4 - y3])

    test1 = np.abs(np.cross(r, s))
    test2 = np.abs(np.cross((q - p), r))

    t0 = np.dot(q - p, r) / np.dot(r, r)
    t1 = t0 + np.dot(s, r) / np.dot(r, r)
    t2 = np.dot(p - q, s) / np.dot(s, s)
    t3 = t2 + np.dot(r, s) / np.dot(s, s)


    inrange = lambda x: x > (0 - tol) and x < (1 + tol)
    distance = lambda x, y: np.sqrt((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2)

    if test1 < tol:

        if test2 < tol:

            if inrange(t0):
                p1 = tuple(p + t0 * r)

            if inrange(t1):
                p2 = tuple(p + t1 * r)

            if inrange(t2):
                p1 = tuple(q + t2 * s)

            if inrange(t3):
                p2 = tuple(q + t3 * s)

    else:
        up = (-x1 * y2 + x1 * y3 + x2 * y1 - x2 * y3 - x3 * y1 + x3 * y2) / (
            x1 * y3 - x1 * y4 - x2 * y3 + x2 * y4 - x3 * y1 + x3 * y2 + x4 * y1 - x4 * y2
        )
        tp = (x1 * y3 - x1 * y4 - x3 * y1 + x3 * y4 + x4 * y1 - x4 * y3) / (
            x1 * y3 - x1 * y4 - x2 * y3 + x2 * y4 - x3 * y1 + x3 * y2 + x4 * y1 - x4 * y2
        )
        if inrange(tp) and inrange(up):
            p1 = tuple(p + tp * r)

    if (p1 is not None) and (p2 is not None) and (distance(p1, p2) < tol):
        p2 = None

    return p1, p2


g = Geometry()
g.import_svg("overlap.svg")
N = len(g.lines)
plt.figure(figsize=(12, 12))
counter = 0
for i in range(N):
    line_1 = g.lines[i]
    plt.plot([line_1.start_pt.x, line_1.end_pt.x], [line_1.start_pt.y, line_1.end_pt.y], 'b-')

    for j in range(i+1, N):
        line_2 = g.lines[j]
        p1, p2 = get_intersetion_points(line_1, line_2)

        if p1 is not None:
            plt.scatter(p1[0], p1[1], c="r", marker="o", s=100)
            counter += 1
            # print(p1)

        if p2 is not None:
            plt.scatter(p2[0], p2[1], c="r", marker="o", s=100)
            counter += 1
            # print(p2)


plt.axis("equal")
#plt.legend()
plt.savefig("overlap.png", dpi=330)

plt.cla()
plt.clf()
plt.close()