from pathlib import Path


def get_next_line():
    with open(Path(__file__).parent / "pareto_front.csv", 'r') as f:
        yield from f

def get_next_solution():
    for line in get_next_line():
        line = line.strip().split(',')
        platform = line.pop(0)
        line = [float(vi) for vi in line]
        f1 = line.pop(0)
        f2 = line.pop(0)
        f3 = line.pop(0)
        nodes = line.pop(0)
        r = line.copy()
        yield tuple(r)



mirrorfactors = []

for sol in get_next_solution():
    left = sol[:10]
    right = sol[10:]
    mirrorfactor = 0.0
    for li, ri in zip(left, right):
        diff = abs(li-ri)
        if diff < 0.5:
            diff = 0.0

        mirrorfactor += diff

    mirrorfactors.append(mirrorfactor)


import matplotlib.pyplot as plt

plt.plot(sorted(mirrorfactors), 'b-')
plt.show()