import re
from math import atan2, degrees, fmod, pi

from digital_twin_distiller.geometry import Geometry
from digital_twin_distiller.objects import CircleArc, Line, Node

# set up regular expressions
# Point(2) = {.1, 0,  0, lc}; -> "2", ".1, 0,  0, lc".split(",")[:2]
# Line(4) = {4, 1}; -> "4", "4, 1".split(",") -> item_data = [int(it.strip()) for it in item_data]
# Circle(60) = {33, 291, 57};

rx_dict = {
    "comment": re.compile(r"[/]{2}"),
    "point": re.compile(r"Point\((\d+)\)[\s=]{1,5}\{(.+?)\};"),
    "line": re.compile(r"Line\((\d+)\)[\s=]{1,5}\{(.+?)\};"),
    "circle": re.compile(r"Circle\((\d+)\)[\s=]{1,5}\{(.+?)\};"),
}


def _parse_line(line):
    """
    Do a regex search against all defined regexes and return the key and match result of the first matching regex,
    except the line is commented out.
    """

    for key, rx in rx_dict.items():
        # find all of the matched solutions in a row and give them back as a list of the elements
        for match in rx.finditer(line):
            if match:
                # Firstly, check if the expression is commented out and do not read the line
                if key == "comment":
                    return None, None
                else:
                    return key, match
    return None, None


def _order_circle_arc_points(start: Node, center: Node, end: Node):
    """
    This function reorders the CircleArcs points in a way that the angle
    difference between the start and end point will be less than 180 degrees.

    :param start: The start point of the CircleArc.
    :param center: The center point of the CircleArc.
    :param end: The end point of the CircleArc.

    """

    # Calculate the angles between the x axis and the points.
    theta0 = atan2(start.y - center.y, start.x - center.x)
    theta1 = atan2(end.y - center.y, end.x - center.x)

    # convert the angle range from [-pi/2, pi/2] -> [0, 360)
    theta0 = degrees(fmod(2 * pi + theta0, 2 * pi))
    theta1 = degrees(fmod(2 * pi + theta1, 2 * pi))

    # calculate the angle difference
    dtheta = theta1 - theta0

    # if dtheta < 0 then the start point has bigger angle therefore the points needs to
    # be swapped.
    if dtheta < 0:
        theta0, theta1 = theta1, theta0
        start, end = end, start

    # Swap is needed again if the angle difference is bigger than 180 degrees.
    if abs(dtheta) > 180:
        theta0, theta1 = theta1, theta0
        start, end = end, start

    return start, center, end


def geo_parser(geo_file):
    """
    Collects the imported entities into a geometry object and returns with a geo object of these geometries.

    :param geo_file: the input data file of the given geometry.
    :return: a Geoemtry object with the collected entities.
    """
    geo = Geometry()

    # open the file and read through it line by line
    with open(geo_file) as file_object:
        row = file_object.readline()
        while row:
            # at each line check for a match with a regex
            parts = row.split(";")

            for row_part in parts:
                row_part += ";"  # inserting back the cutted ;

                key, match = _parse_line(row_part)
                if key == "point":
                    # Example:
                    # Point(2) = {.1, 0,  0, lc};
                    item_id = match.group(1)  # "2"
                    item_data = match.group(2)  # ".1, 0,  0, lc"
                    item_data = item_data.split(",")[:2]
                    item_data = [float(it.strip()) for it in item_data]  # [0.1, 0.0]

                    try:
                        geo.add_node(Node(id_=int(item_id), x=float(item_data[0]), y=float(item_data[1])))
                    except:
                        raise ValueError("String data inserted instead of the point data")

                if key == "line":
                    # Example
                    # Line(4) = {4, 1};
                    item_id = match.group(1)  # "4"
                    item_data = match.group(2)  # "4, 1"
                    item_data = item_data.split(",")[:2]
                    item_data = [float(it.strip()) for it in item_data]  # [0.1, 0.0]

                    try:
                        geo.add_line(
                            Line(
                                id_=int(item_id),
                                start_pt=geo.find_node(int(item_data[0])),
                                end_pt=geo.find_node(int(item_data[1])),
                            )
                        )
                    except:
                        raise ValueError("Invalid point objects in the geometry")

                if key == "circle":
                    # Example
                    #  Circle(60) = {33, 291, 57};
                    item_id = match.group(1)  # "60"
                    item_data = match.group(2)  # "33, 291, 57"
                    item_data = item_data.split(",")[:3]
                    item_data = [float(it.strip()) for it in item_data]  # [0.1, 0.0]

                    try:
                        start, center, end = _order_circle_arc_points(
                            geo.find_node(int(item_data[0])),
                            geo.find_node(int(item_data[1])),
                            geo.find_node(int(item_data[2])),
                        )
                        geo.add_arc(CircleArc(id_=int(item_id), start_pt=start, center_pt=center, end_pt=end))
                    except:
                        raise ValueError("Invalid point objects in the geometry")

            row = file_object.readline()

    return geo
