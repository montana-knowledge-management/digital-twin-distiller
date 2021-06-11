# TODO: remove matplotlib after testing
import numpy as np
import matplotlib.pyplot as plt
from adze_modeler.geometry import Geometry
g = Geometry()
g.import_svg("overlap.svg")
g.sanitize_geometry()
plt.figure(figsize=(12, 12))
for i, li in enumerate(g.lines):
    plt.plot([li.start_pt.x, li.end_pt.x], [li.start_pt.y, li.end_pt.y])


plt.axis("equal")
# plt.legend()

plt.savefig(f"temp.png", dpi=330)
plt.cla()
plt.clf()
plt.close()