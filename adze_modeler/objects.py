import math


class Node:
    """
    A Node identified by (x,y) coordinates, optionally it can contains an id number or a label. The id_number and
    the label can be important to rotate and copy and rotate the selected part of the geometry.
    """

    def __init__(self, x=0.0, y=0.0, id=None, label=None, precision=6):
        self.x = x
        self.y = y
        self.id = id  # a node has to got a unique id to be translated or moved
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
        return f"({self.x}, {self.y}, id={self.id},label={self.label})"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.x!r}, {self.y!r}, id={self.id!r},label={self.label!r})"

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
        self.id = id
        self.label = label

    def clone(self):
        return Line(self.start_pt, self.end_pt, self.id, self.label)

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
        return f"{self.__class__.__name__}({self.start_pt!r}, {self.end_pt!r}, id={self.id!r},label={self.label!r})"


class CircleArc:
    """A directed line, which is defined by the (start -> end) points"""

    def __init__(self, start_pt, center_pt, end_pt, id=None, label=None):

        self.start_pt = start_pt
        self.center_pt = center_pt
        self.end_pt = end_pt
        self.id = id
        self.label = label

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
