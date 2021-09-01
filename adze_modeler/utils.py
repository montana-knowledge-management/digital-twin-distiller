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


if __name__ == "__main__":
    from adze_modeler.objects import Node
    import matplotlib.pyplot as plt

    p1 = Node(2, 2)
    p2 = Node(5, 5)
    p3 = Node(-3, 3)

    p4 = mirror_point(p1, p2, p3)
    print(p4)

    plt.plot([p1.x, p2.x], [p1.y, p2.y], "k")
    plt.scatter(*p3, c="b")
    plt.scatter(*p4, c="b")
    plt.axis("equal")
    plt.show()
