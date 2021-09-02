from uuid import uuid4


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
