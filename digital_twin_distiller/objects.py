import math
from collections.abc import Iterable
from copy import copy

from numpy import linspace

from digital_twin_distiller.utils import getID, get_phi, get_short_id, mirror_point, pairwise


def transformIntoInterval(minx, maxx, x):
    return min(max(x, minx), maxx)


def transformMeshScalingInterval(x):
    return transformIntoInterval(0.1, 20, x)


def calculateMeshScalingInput(meshScaling):
    return abs(int(math.ceil(float(meshScaling))))


# def transform(x):
#     a = 0.1
#     b = 20
#     if a <= x <= b:
#         return x
#     # [0.1, num] => [0.1, 20]
#     else:
#         res = (x - 0.1) * ((b - a) / x) + a
#         if res > a:
#             return res
#         return a


class Node:
    """
    A Node identified by (x,y) coordinates, optionally it can contains an id number or a label. The id_number and
    the label can be important to rotate and copy and rotate the selected part of the geometry.
    """

    def __init__(self, x=0.0, y=0.0, id_=None, label=None, precision=6):
        self.x = x
        self.y = y
        self.id = id_ or getID()  # a node has to got a unique id to be translated or moved
        self.label = label  # can be used to denote a group of the elements and make some operation with them
        self.precision = precision  # number of the digits, every coordinate represented in the same precision
        self.hanging = True  # if its contained by another object it will be set to False

    @classmethod
    def from_polar(cls, r: float, phi: float):
        """
        Create a Node from its polar coordinates.

        :param float r: the length of the vector
        :param float phi: the angle of the vector in degrees
        """
        return cls(r * math.cos(math.radians(phi)), r * math.sin(math.radians(phi)))

    def __getitem__(self, item):
        if item == 0:
            return self.x
        elif item == 1:
            return self.y
        else:
            raise IndexError

    def __le__(self, other):
        if self < other or self == other:
            return True

        return False

    def __ge__(self, other):
        if self > other or self == other:
            return True

        return False

    def __gt__(self, other):
        if not self == other and not self < other:
            return True
        else:
            return False

    def __eq__(self, other):
        return abs(self.x - other.x) < 1e-5 and abs(self.y - other.y) < 1e-5

    def __add__(self, p):
        """Point(x1+x2, y1+y2)"""
        if isinstance(p, Node):
            return Node(self.x + p.x, self.y + p.y)
        elif isinstance(p, Iterable):
            return Node(*(pi + ci for pi, ci in zip(p, self)))
        else:
            return Node(self.x + p, self.y + p)

    def __sub__(self, p):
        """Point(x1-x2, y1-y2)"""
        return Node(self.x - p.x, self.y - p.y)

    def __mul__(self, scalar):
        """Point(x1*x2, y1*y2)"""
        return Node(self.x * scalar, self.y * scalar)

    def __rmul__(self, scalar):
        return self * scalar

    def __truediv__(self, scalar):
        return Node(self.x / scalar, self.y / scalar)

    def __matmul__(self, other):
        """
        Dot prduct
        n1 @ n2
        """
        return self.x * other.x + self.y * other.y

    def __str__(self):
        # return f"({self.x:.1f}, {self.y:.1f}, label={self.label})"
        return f"{self.__class__.__name__}({self.x:.1f}, {self.y:.1f}, id={hex(self.id)[-5:]})"
        # return f"N({self.x:.1f}, {self.y:.1f}, {self.label})"
        # return f"{self.label}"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.x!r}, {self.y!r}, id={self.id!r},label={self.label!r})"
        # return f"N({self.x:.1f}, {self.y:.1f}, {self.label})"
        # return f"{self.label}"

    def __copy__(self):
        return Node(
            self.x,
            self.y,
            id_=getID(),
            label=self.label,
            precision=self.precision,
        )

    def __iter__(self):
        yield from (self.x, self.y)

    def __abs__(self):
        return self.length()

    def __lt__(self, o):
        """
        This function compares the operands l2 value with its l2 value. If they're equal (whitin tolerance)
        that means the nodes are on the same circle. If that's the case then check their angle return accordingly.
        If the l2 values are different then one node is obviously bigger than the other. In that case check the
        differences sign.
        :param o: other Node
        :return: True or False
        """
        diff = abs(self) - abs(o)
        tol = 1e-5

        phi_self = get_phi(*self, tol=tol)
        phi_other = get_phi(*o, tol=tol)
        angle_diff = phi_self - phi_other

        if abs(angle_diff) < tol:
            if diff < -tol:
                return True
            else:
                return False
        elif angle_diff < -tol:
            return True
        else:
            return False

        # if abs(diff) < tol:
        #     # They're on the same circle, more checks needed to decide
        #     if (phi_self - phi_other) < - tol:
        #         return True
        #     else:
        #         return False
        # elif diff < -tol:
        #     return True
        # else:
        #     return False

    def __hash__(self):
        return int(get_short_id(self), base=16)

    def length(self):
        return math.hypot(*self)

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
        s, c = (f(rad) for f in (math.sin, math.cos))
        x, y = (c * self.x - s * self.y, s * self.x + c * self.y)
        return Node(round(x, self.precision), round(y, self.precision))

    def rotate_about(self, p, theta):
        """Rotate counter-clockwise around a point, by theta radians. The new position is returned as a new Point."""
        result = self.clone()
        result.move_xy(-p.x, -p.y)
        result = result.rotate(theta)
        result.move_xy(p.x, p.y)
        return result

    def unit_to(self, other):
        """
        This function returns a unit vector that points from self to other
        """
        u = other - self
        return u * (1 / u.length())

    def angle_to(self, other):
        """
        This function returns the angle between self and an another Node
        instance.
        """
        return math.acos((self @ other) / (self.length() * other.length()))

    def mean(self, other):
        return (self + other) / 2


class Line:
    """A directed line, which is defined by the (start -> end) points"""

    def __init__(self, start_pt, end_pt, id_=None, label=None, color=None, attributes: dict = {}):
        # sorting the incoming points by coordinate
        # sorted_points = sorted((start_pt, end_pt), key=lambda pi: pi.x)  # sorting by x coordinate
        # sorted_points = sorted(sorted_points, key=lambda pi: pi.y)  # sorting by y coordinate
        # self.start_pt = sorted_points[0]
        # self.end_pt = sorted_points[-1]
        self.start_pt = start_pt
        self.end_pt = end_pt
        self.id = id_ or getID()
        self.label = label
        self.color = color  # the color of the given edge can be used to render the appropriate boundary conditions to the given edges
        self.attributes = attributes.copy()
        self.length = math.dist(self.start_pt, self.end_pt)
        self.meshScaling = 1.0

        getMeshScaling = attributes.get("meshScaling")

        if getMeshScaling:
            self.set_mesh_scaling(getMeshScaling)

    def __copy__(self):
        return Line(
            copy(self.start_pt),
            copy(self.end_pt),
            id_=getID(),
            label=self.label,
            color=self.color,
            attributes=self.attributes,
        )

    def set_mesh_scaling(self, meshScaling):
        self.meshScaling = transformMeshScalingInterval(self.length / calculateMeshScalingInput(meshScaling))

    def distance_to_point(self, px, py):
        """
        This function calculates the minimum distance between a line segment and a point
        https://www.geeksforgeeks.org/minimum-distance-from-a-point-to-the-line-segment-using-vectors/
        """
        # p = Node(px, py)
        # center_pt = (self.start_pt + self.end_pt) / 2
        # d1 = self.start_pt.distance_to(p)
        # d2 = center_pt.distance_to(p)
        # d3 = self.end_pt.distance_to(p)
        # return min(d1, d2, d3)

        A = (self.start_pt.x, self.start_pt.y)
        B = (self.end_pt.x, self.end_pt.y)
        E = (px, py)

        # vector AB
        AB = [None, None]
        AB[0] = B[0] - A[0]
        AB[1] = B[1] - A[1]

        # vector BP
        BE = [None, None]
        BE[0] = E[0] - B[0]
        BE[1] = E[1] - B[1]

        # vector AP
        AE = [None, None]
        AE[0] = E[0] - A[0]
        AE[1] = E[1] - A[1]

        # Variables to store dot product

        # Calculating the dot product
        AB_BE = AB[0] * BE[0] + AB[1] * BE[1]
        AB_AE = AB[0] * AE[0] + AB[1] * AE[1]

        # Minimum distance from
        # point E to the line segment
        reqAns = 0

        # Case 1
        if AB_BE > 0:

            # Finding the magnitude
            y = E[1] - B[1]
            x = E[0] - B[0]
            reqAns = math.sqrt(x * x + y * y)

        # Case 2
        elif AB_AE < 0:
            y = E[1] - A[1]
            x = E[0] - A[0]
            reqAns = math.sqrt(x * x + y * y)

        # Case 3
        else:

            # Finding the perpendicular distance
            x1 = AB[0]
            y1 = AB[1]
            x2 = AE[0]
            y2 = AE[1]
            mod = math.sqrt(x1 * x1 + y1 * y1)
            reqAns = abs(x1 * y2 - y1 * x2) / mod

        return reqAns

    def __eq__(self, other):
        """
        TODO: docstring here
        """
        d1 = self.start_pt.distance_to(other.start_pt)
        d2 = self.start_pt.distance_to(other.end_pt)
        d3 = self.end_pt.distance_to(other.start_pt)
        d4 = self.end_pt.distance_to(other.end_pt)
        distances = sorted([d1, d2, d3, d4])
        if distances[0] < 1e-5 and distances[1] < 1e-5:
            return True
        else:
            return False

    def __call__(self, t: float):
        assert (0 <= t) and (t <= 1), f"t [0, 1] not {t}"
        return self.start_pt + (self.end_pt - self.start_pt) * t

    def __repr__(self):
        return f"{self.__class__.__name__}({self.start_pt}, {self.end_pt},label={self.label!r}, color={self.color})"
        # return f"{self.__class__.__name__}({self.start_pt}, {self.end_pt}, id={hex(self.id)[-5:]})"


class CircleArc:
    """A directed line, which is defined by the (start -> end) points"""

    def __init__(
        self, start_pt, center_pt, end_pt, id_=None, label=None, max_seg_deg=1, color=None, attributes: dict = {}
    ):
        self.start_pt = start_pt
        self.center_pt = center_pt
        self.end_pt = end_pt
        self.id = id_ or getID()
        self.label = label
        self.max_seg_deg = max_seg_deg
        self.color = color
        self.attributes = attributes.copy()

        self.meshScaling = attributes.get("meshScaling")

        self.radius = self.start_pt.distance_to(self.center_pt)
        clamp = self.start_pt.distance_to(self.end_pt) / 2.0
        try:
            self.theta = round(math.asin(clamp / self.radius) * 180 / math.pi * 2, 2)
            self.apex_pt = self.start_pt.rotate_about(self.center_pt, math.radians(self.theta / 2))

            if self.meshScaling:
                self.set_max_seg_deg_by_mesh_scaling(self.meshScaling)

        except ValueError:
            self.apex_pt = self.start_pt.rotate_about(self.center_pt, math.radians(90))
            self.theta = 180

    def set_max_seg_deg_by_mesh_scaling(self, meshScaling):
        self.max_seg_deg = transformMeshScalingInterval(self.theta / calculateMeshScalingInput(meshScaling))

    @classmethod
    def from_radius(cls, start_pt: Node, end_pt: Node, r: float = 1.0, attributes: dict = {}):
        """
        Construct a CircleArc instance from start- and end-point and a radius.

        Parameters:
            start_pt: Starting point of the arc.
            end_pt: End point of the arc.
            r: radius of the arc. r > 0.
            attributes: dict
        """
        assert r > 0.0, f"Radius should be greater than 0. Got {r=}."
        assert start_pt != end_pt, "Start- and End-point cannot be the same."

        a = start_pt.distance_to(end_pt)
        alpha = math.acos(a * 0.5 / r)
        u = start_pt.unit_to(end_pt) * r
        u = u.rotate(alpha) + start_pt
        return cls(start_pt, u, end_pt, attributes=attributes)

    def distance_to_point(self, x, y):
        """
        This function returns the minimum distance between p and the circle arcs points:
        start, end, center and apex point.
        """
        p = Node(x, y)
        d1 = self.start_pt.distance_to(p)
        d2 = self.apex_pt.distance_to(p)
        d3 = self.end_pt.distance_to(p)
        return min(d1, d2, d3)

    def __eq__(self, other):
        """
        If 2 circles have the same set of points, then they are equal. Any difference will result a False
        return value.
        """

        if self.start_pt != other.start_pt:
            return False

        if self.center_pt != other.center_pt:
            return False

        if self.apex_pt != other.apex_pt:
            return False

        if self.end_pt != other.end_pt:
            return False

        return True

    def __copy__(self):
        return CircleArc(
            copy(self.start_pt),
            copy(self.center_pt),
            copy(self.end_pt),
            max_seg_deg=self.max_seg_deg,
            color=self.color,
            attributes=self.attributes,
        )

    def __repr__(self):
        return "{}({!r}, {!r}, {!r}, id={!r},label={!r}, color={!r})".format(
            self.__class__.__name__, self.start_pt, self.center_pt, self.end_pt, self.id, self.label, self.color
        )


class CubicBezier:
    def __init__(
        self,
        start_pt,
        control1,
        control2,
        end_pt,
        id_=None,
        label=None,
        color=None,
        attributes: dict = {},
        n_segment=51,
    ):
        self.start_pt = start_pt
        self.control1 = control1
        self.control2 = control2
        self.end_pt = end_pt
        self.id = id_ or getID()
        self.label = label
        self.color = color
        self.attributes = attributes.copy()
        self.n_segment = attributes.get("meshScaling", n_segment)

    def approximate(self):
        X, Y = zip(*(self(ti) for ti in linspace(0, 1, self.n_segment + 1)))

        for Xi, Yi in zip(pairwise(X), pairwise(Y)):
            n0 = Node(Xi[0], Yi[0])
            n1 = Node(Xi[1], Yi[1])
            yield Line(n0, n1)

    def __call__(self, t: float):
        assert (0 <= t) and (t <= 1), f"t [0, 1] not {t}"
        X = (
            (1 - t) ** 3 * self.start_pt.x
            + 3 * (1 - t) ** 2 * t * self.control1.x
            + 3 * (1 - t) * t**2 * self.control2.x
            + t**3 * self.end_pt.x
        )

        Y = (
            (1 - t) ** 3 * self.start_pt.y
            + 3 * (1 - t) ** 2 * t * self.control1.y
            + 3 * (1 - t) * t**2 * self.control2.y
            + t**3 * self.end_pt.y
        )

        return X, Y

    def __eq__(self, other):
        """
        If 2 Bezier-Curves have the same set of points, then they are equal.
        """

        if math.dist(self.start_pt, other.start_pt) > 1e-5:
            return False

        if math.dist(self.control1, other.control1) > 1e-5:
            return False

        if math.dist(self.control2, other.control2) > 1e-5:
            return False

        if math.dist(self.end_pt, other.end_pt) > 1e-5:
            return False

        return True

    def __repr__(self):
        return "{}({!r}, {!r}, {!r}, {!r}, id={!r},label={!r}, color={!r})".format(
            self.__class__.__name__,
            self.start_pt,
            self.control1,
            self.control2,
            self.end_pt,
            self.id,
            self.label,
            self.color,
        )


class Rectangle:
    def __init__(self, x0: float = 0.0, y0: float = 0.0, **kwargs):
        """

         d --------------------------- [c]
         |                             |
         |                             |
         |                             |
        [a] -------------------------- b

        a: [x0, y0] is the fixed point of the rectangle
        c: [x1, y1] is the upper right corner of the rectangle

        keyword arguments:
            - width, height: width=length of the a-b/d-c line, height=length of the a-d/b-c line, width and height can
            be negative
            - x1, y1: specifies the c point

        """
        # coordinates of the 4 points
        self.a = Node(x0, y0)
        self.b = None
        self.c = None
        self.d = None

        # center-point of the rectangle
        self.cp = None

        # length of the line a-b and d-c
        self.width = 0.0

        # length of the line a-d and b-c
        self.height = 0.0

        if {"width", "height"}.issubset(kwargs.keys()):
            w = kwargs["width"]
            h = kwargs["height"]

            self.b = Node(self.a.x + w, self.a.y)
            self.c = Node(self.a.x + w, self.a.y + h)
            self.d = Node(self.a.x, self.a.y + h)

        elif {"x1", "y1"}.issubset(kwargs.keys()):
            self.c = Node(kwargs["x1"], kwargs["y1"])
            self.b = Node(self.c.x, self.a.y)
            self.d = Node(self.a.x, self.c.y)

        else:
            raise ValueError("Not enough parameters were given.")

        self._calc_centerpoint()
        self._calc_width_height()

    def rotate(self, phi: float, fx_point=None):
        """
        Rotate a Rectangle instance around a point with phi degrees. The default point is the center-point.

        :param fx_point: Sets one of the points as the origin of the rotation. Accepted values are 'a', 'b', 'c', 'd'
        """
        phi = phi * math.pi / 180
        rotation_center = self.cp

        if fx_point is not None:
            if fx_point == "a":
                rotation_center = self.a
            elif fx_point == "b":
                rotation_center = self.b
            elif fx_point == "c":
                rotation_center = self.c
            elif fx_point == "d":
                rotation_center = self.d
            else:
                raise ValueError(f"Invalid value for fx_point. Got {fx_point=}")

        self.a = self.a.rotate_about(rotation_center, phi)
        self.b = self.b.rotate_about(rotation_center, phi)
        self.c = self.c.rotate_about(rotation_center, phi)
        self.d = self.d.rotate_about(rotation_center, phi)

        self._calc_centerpoint()

    def set_width(self, new_width, fx_point=None):
        """
        Sets the width of the rectangle from a fixed point. This point is the center-point by default.
        """
        difference = new_width - self.width

        # unit vectors
        u_ab = self.a.unit_to(self.b)
        u_dc = self.d.unit_to(self.c)

        if fx_point is None:
            self.a = self.a - u_ab * difference / 2
            self.b = self.b + u_ab * difference / 2
            self.d = self.d - u_dc * difference / 2
            self.c = self.c + u_dc * difference / 2
        elif fx_point in {"a", "d"}:
            self.b = self.b + u_ab * difference
            self.c = self.c + u_dc * difference
        elif fx_point in {"c", "b"}:
            self.a = self.a - u_ab * difference
            self.d = self.d - u_dc * difference
        else:
            raise ValueError(f"Invalid value for fx_point. Got {fx_point=}")

        self._calc_width_height()

    def set_height(self, new_height, fx_point=None):
        """
        Sets the height of the rectangle from a fixed point. This point is the center-point by default.
        Other points of the rectangle can be used as a reference.

        :param new_height: the new height of the rectangle
        :param fx_point: 'a', 'b', 'c', 'd'
        """
        difference = new_height - self.height
        # unit vectors
        u_ad = self.a.unit_to(self.d)
        u_bc = self.b.unit_to(self.c)

        if fx_point is None:
            self.a = self.a - u_ad * difference / 2
            self.b = self.b - u_ad * difference / 2
            self.c = self.c + u_bc * difference / 2
            self.d = self.d + u_bc * difference / 2
        elif fx_point in {"a", "b"}:
            self.c = self.c + u_bc * difference
            self.d = self.d + u_ad * difference
        elif fx_point in {"c", "d"}:
            self.a = self.a - u_ad * difference
            self.b = self.b - u_bc * difference
        else:
            raise ValueError(f"Invalid value for fx_point. Got {fx_point=}")

        self._calc_width_height()

    def put(self, x, y, fx_point=None):
        """
        This function moves the rectangle such that the fx_point touches the (x, y) point. The default fx_point is
        the center-point.
        """
        ref = Node(x, y)
        difference = None

        if fx_point is None:
            difference = ref - self.cp
        elif fx_point == "a":
            difference = ref - self.a
        elif fx_point == "b":
            difference = ref - self.b
        elif fx_point == "c":
            difference = ref - self.c
        elif fx_point == "d":
            difference = ref - self.d
        else:
            raise ValueError(f"Invalid value for fx_point. Got {fx_point=}")

        self.translate(difference.x, difference.y)

    def translate(self, dx=0, dy=0):
        self.a.move_xy(dx, dy)
        self.b.move_xy(dx, dy)
        self.c.move_xy(dx, dy)
        self.d.move_xy(dx, dy)
        self.cp.move_xy(dx, dy)

    def mirror(self, p1=(0, 0), p2=(0, 1)):
        p1 = Node(*p1)
        p2 = Node(*p2)
        self.a = mirror_point(p1, p2, self.a)
        self.b = mirror_point(p1, p2, self.b)
        self.c = mirror_point(p1, p2, self.c)
        self.d = mirror_point(p1, p2, self.d)
        self._calc_centerpoint()

    def _calc_centerpoint(self):
        """Calculates the center-point of the rectangle."""
        self.cp = Node(
            (self.a.x + self.b.x + self.c.x + self.d.x) / 4,
            (self.a.y + self.b.y + self.c.y + self.d.y) / 4,
        )

    def _calc_width_height(self):
        """Calculates the width and the height of the Rectangle."""
        self.width = math.hypot(self.b.x - self.a.x, self.b.y - self.a.y)
        self.height = math.hypot(self.d.x - self.a.x, self.d.y - self.a.y)
        self._calc_centerpoint()

    def _print_sidelengths(self):
        """
        This function is for debugging purposes only. It prints out the side lengths of the rectangle.
        """

        print("--" * 15)
        print("ab:", math.hypot(self.b.x - self.a.x, self.b.y - self.a.y))
        print("dc:", math.hypot(self.d.x - self.c.x, self.d.y - self.c.y))
        print()
        print("ad:", math.hypot(self.d.x - self.a.x, self.d.y - self.a.y))
        print("bc:", math.hypot(self.c.x - self.b.x, self.c.y - self.b.y))
        print("--" * 15)

    def __iter__(self):
        yield from (
            self.a.x,
            self.a.y,
            self.b.x,
            self.b.y,
            self.c.x,
            self.c.y,
            self.d.x,
            self.d.y,
        )

    def __copy__(self):
        r = Rectangle(self.a.x, self.a.y, width=1.0, height=1.0)
        r.a = copy(self.a)
        r.b = copy(self.b)
        r.c = copy(self.c)
        r.d = copy(self.d)
        r._calc_width_height()
        r._calc_centerpoint()
        return r

    def __repr__(self):
        return (
            f"a:({self.a.x:.2f}, {self.a.y:.2f}) cp:({self.cp.x:.2f}, {self.cp.y:.2f}) "
            f"c: ({self.c.x:.2f}, {self.c.y:.2f})"
            f" w={self.width:.3f}, height={self.height:.3f}"
        )
