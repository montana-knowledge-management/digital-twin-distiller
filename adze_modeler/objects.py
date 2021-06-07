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
        self.start_pt = start_pt
        self.end_pt = end_pt
        self.id = id
        self.label = label

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
