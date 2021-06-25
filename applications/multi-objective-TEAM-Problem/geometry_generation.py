import random

import svgpathtools as svg


def get_rectangle(h, w, radii, offsetz, printlabel=True):
    l1 = svg.Line(0 + 0j, w + 0j)
    l2 = svg.Line(w + 0j, w + h * 1j)
    l3 = svg.Line(w + h * 1j, 0 + h * 1j)
    l4 = svg.Line(0 + h * 1j, 0 + 0j)

    rectangle = svg.Path(l1, l2, l3, l4)
    rectangle = rectangle.translated(radii + offsetz * 1j)

    xmin, xmax, ymin, ymax = rectangle.bbox()
    label_x = (xmax - xmin) / 2 + xmin
    label_y = (ymax - ymin) / 2 + ymin
    if printlabel:
        print(f"snapshot.assign_material({label_x:.3f}, {label_y:.3f}, name='J+')")

    return rectangle


random.seed(42)

# Number of coils
N = 20

# height of the coil
h = 1.5

# width of the coil
w = 1

paths = []
for i in range(N):
    R = random.uniform(5, 50)
    offsetz = i * h - N * h / 2
    paths.append(get_rectangle(h, w, R, offsetz))


paths.append(get_rectangle(10, 5, 0, -5, printlabel=False))
paths.append(get_rectangle(70, 70, 0, -35, printlabel=False))
svg.wsvg(paths, filename="geometry.svg")
