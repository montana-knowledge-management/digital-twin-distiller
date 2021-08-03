import math
from collections.abc import Iterable
from copy import copy

from adze_modeler.utils import getID


class Node:
    """
    A Node identified by (x,y) coordinates, optionally it can contains an id number or a label. The id_number and
    the label can be important to rotate and copy and rotate the selected part of the geometry.
    """

    def __init__(self, x=0.0, y=0.0, id=None, label=None, precision=6):
        self.x = x
        self.y = y
        self.id = id or getID()  # a node has to got a unique id to be translated or moved
        self.label = label  # can be used to denote a group of the elements and make some operation with them
        self.precision = precision  # number of the digits, every coordinate represented in the same precision
        self.hanging = True  # if its contained by another object it will be set to False

    def __add__(self, p):
        """Point(x1+x2, y1+y2)"""
        return Node(self.x + p.x, self.y + p.y)

    def __sub__(self, p):
        """Point(x1-x2, y1-y2)"""
        return Node(self.x - p.x, self.y - p.y)

    def __mul__(self, scalar):
        """Point(x1*x2, y1*y2)"""
        return Node(self.x * scalar, self.y * scalar)

    def __rmul__(self, other):
        """
        Dot prduct
        n1 @ n2
        """
        return self.x * other.x + self.y * other.y

    def __str__(self):
        return f"({self.x:.1f}, {self.y:.1f}, label={self.label})"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.x!r}, {self.y!r}, id={self.id!r},label={self.label!r})"

    def __copy__(self):
        return Node(self.x, self.y, id=getID(), label=self.label, precision=self.precision)

    def length(self):
        return math.sqrt(self.x ** 2 + self.y ** 2)

    def distance_to(self, p):
        """Calculate the distance between two points."""
        return (self - p).length()

    def as_tuple(self):
        """(x, y)"""
        return (self.x, self.y)

    def clone(self):
        """Return a full copy of this point."""
        return Node(self.x, self.y, self.id, self.label, self.precision)

    def move_xy(self, dx, dy):
        """Move to new (x+dx,y+dy)."""
        self.x = round(self.x + dx, self.precision)
        self.y = round(self.y + dy, self.precision)

    def rotate(self, rad):
        """Rotate counter-clockwise by rad radians.

        Positive y goes *up,* as in traditional mathematics.

        Interestingly, you can use this in y-down computer graphics, if
        you just remember that it turns clockwise, rather than
        counter-clockwise.

        The new position is returned as a new Point.
        """
        s, c = [f(rad) for f in (math.sin, math.cos)]
        x, y = (c * self.x - s * self.y, s * self.x + c * self.y)
        return Node(round(x, self.precision), round(y, self.precision))

    def rotate_about(self, p, theta):
        """Rotate counter-clockwise around a point, by theta degrees. The new position is returned as a new Point."""
        result = self.clone()
        result.move_xy(-p.x, -p.y)
        result = result.rotate(theta)
        result.move_xy(p.x, p.y)
        return result


class Line:
    """A directed line, which is defined by the (start -> end) points"""

    def __init__(self, start_pt, end_pt, id=None, label=None):
        # sorting the incoming points by coordinate
        # sorted_points = sorted((start_pt, end_pt), key=lambda pi: pi.x)  # sorting by x coordinate
        # sorted_points = sorted(sorted_points, key=lambda pi: pi.y)  # sorting by y coordinate
        # self.start_pt = sorted_points[0]
        # self.end_pt = sorted_points[-1]
        self.start_pt = start_pt
        self.end_pt = end_pt
        self.id = id or getID()
        self.label = label

    def __copy__(self):
        return Line(copy(self.start_pt), copy(self.end_pt), id=getID(), label=self.label)

    def distance_to_point(self, px, py):
        """
        This function calculates the minimum distance between a line segment and a point
        https://www.geeksforgeeks.org/minimum-distance-from-a-point-to-the-line-segment-using-vectors/
        """
        A = (self.start_pt.x, self.start_pt.y)
        B = (self.end_pt.x, self.end_pt.y)
        E = (px, py)

        # vector AB
        AB = [None, None];
        AB[0] = B[0] - A[0];
        AB[1] = B[1] - A[1];

        # vector BP
        BE = [None, None];
        BE[0] = E[0] - B[0];
        BE[1] = E[1] - B[1];

        # vector AP
        AE = [None, None];
        AE[0] = E[0] - A[0];
        AE[1] = E[1] - A[1];

        # Variables to store dot product

        # Calculating the dot product
        AB_BE = AB[0] * BE[0] + AB[1] * BE[1];
        AB_AE = AB[0] * AE[0] + AB[1] * AE[1];

        # Minimum distance from
        # point E to the line segment
        reqAns = 0;

        # Case 1
        if (AB_BE > 0):

            # Finding the magnitude
            y = E[1] - B[1];
            x = E[0] - B[0];
            reqAns = math.sqrt(x * x + y * y);

        # Case 2
        elif (AB_AE < 0):
            y = E[1] - A[1];
            x = E[0] - A[0];
            reqAns = math.sqrt(x * x + y * y);

        # Case 3
        else:

            # Finding the perpendicular distance
            x1 = AB[0];
            y1 = AB[1];
            x2 = AE[0];
            y2 = AE[1];
            mod = math.sqrt(x1 * x1 + y1 * y1);
            reqAns = abs(x1 * y2 - y1 * x2) / mod

        return reqAns

    def __repr__(self):
        return f"{self.__class__.__name__}({self.start_pt}, {self.end_pt},label={self.label!r})"


class CircleArc:
    """A directed line, which is defined by the (start -> end) points"""

    def __init__(self, start_pt, center_pt, end_pt, id=None, label=None):
        self.start_pt = start_pt
        self.center_pt = center_pt
        self.end_pt = end_pt
        self.id = id or getID()
        self.label = label
        self.max_seg_deg = 20

    def __copy__(self):
        return CircleArc(copy(self.start_pt), copy(self.center_pt), copy(self.end_pt))

    def __repr__(self):
        return "{}({!r}, {!r}, {!r}, id={!r},label={!r})".format(
            self.__class__.__name__, self.start_pt, self.center_pt, self.end_pt, self.id, self.label
        )


class CubicBezier:
    def __init__(self, start_pt, control1, control2, end_pt, id=None, label=None):
        self.start_pt = start_pt
        self.control1 = control1
        self.control2 = control2
        self.end_pt = end_pt
        self.id = id
        self.label = label

    def __repr__(self):
        return "{}({!r}, {!r}, {!r}, {!r}, id={!r},label={!r})".format(
            self.__class__.__name__, self.start_pt, self.control1, self.control2, self.end_pt, self.id, self.label
        )


class ParametricBezier:
    def __init__(self, **kwargs):
        self.p0 = kwargs.get('start_pt', [None, None])
        self.p1 = kwargs.get('c1', [None, None])
        self.p2 = kwargs.get('c2', [None, None])
        self.p3 = kwargs.get('end_pt', [None, None])

    def set(self, **kwargs):
        self.p0 = kwargs.get('start_pt', self.p0)
        self.p1 = kwargs.get('c1', self.p1)
        self.p2 = kwargs.get('c2', self.p2)
        self.p3 = kwargs.get('end_pt', self.p3)

    def casteljau(self, p0, p1, p2, p3):
        m = ((p1[0] +p2[0])*0.5,
             (p1[1] + p2[1]) * 0.5)
        l0 = p0
        r3 = p3

        l1=((p0[0]+p1[0])*0.5,
            (p0[1]+p1[1])*0.5)
        r2 = ((p2[0]+p3[0])*0.5,
              (p2[1]+p3[1])*0.5)

        l2 = ((l1[0] + m[0]) * 0.5,
              (l1[1] + m[1]) * 0.5)
        r1 = ((r2[0] + m[0]) * 0.5,
              (r2[1] + m[1]) * 0.5)

        l3 = ((l2[0]+r1[0])*0.5,
              (l2[1]+r1[1])*0.5)
        r0 = l3

        return (r0,r1,r2,r3), (l0,l1,l2,l3)

    def approximate(self, nb_iter=0):
        lines = [(self.p0,self.p1, self.p2, self.p3)]
        for iter_i in range(nb_iter):
            templines = []
            for curve_i in lines:
                r, l = self.casteljau(*curve_i)
                templines.append(l)
                templines.append(r)
            lines.clear()
            lines = templines.copy()

        linex = [(ci[0][0], ci[-1][0]) for ci in lines]
        linex = [item for sublist in linex for item in sublist]
        liney = [(ci[0][1], ci[-1][1]) for ci in lines]
        liney = [item for sublist in liney for item in sublist]

        return linex, liney

    def __call__(self, t: float):
        assert (0 <= t) and (t <= 1), f"t [0, 1] not {t}"
        X = (1 - t) ** 3 * self.p0[0] + 3 * (1 - t) ** 2 * t * self.p1[0] + 3 * (1 - t) * t ** 2 * self.p2[0] + t ** 3 * \
            self.p3[0]

        Y = (1 - t) ** 3 * self.p0[1] + 3 * (1 - t) ** 2 * t * self.p1[1] + 3 * (1 - t) * t ** 2 * self.p2[1] + t ** 3 * \
            self.p3[1]

        return X, Y