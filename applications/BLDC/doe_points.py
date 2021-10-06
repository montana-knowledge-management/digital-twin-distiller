import matplotlib.pyplot as plt
import itertools
from model import *
from math import hypot
from adze_modeler.doe import *

def add_scatter_point(ax, x, y, z, color):
    ax.scatter(x, y, z, c=color, s=70, edgecolor='k', alpha=1)

def get_unit_cube():
    # fig = plt.figure()
    ax = plt.axes(projection='3d')
    # plane @ z=-1
    ax.plot3D([-1, 1], [-1, -1], [-1, -1], 'k--')
    ax.plot3D([-1, 1], [1, 1], [-1, -1], 'k--')
    ax.plot3D([-1, -1], [-1, 1], [-1, -1], 'k--')
    ax.plot3D([1, 1], [-1, 1], [-1, -1], 'k--')

    # plane @ z=1
    ax.plot3D([-1, 1], [-1, -1], [1, 1], 'k--')
    ax.plot3D([-1, 1], [1, 1], [1, 1], 'k--')
    ax.plot3D([-1, -1], [-1, 1], [1, 1], 'k--')
    ax.plot3D([1, 1], [-1, 1], [1, 1], 'k--')

    # pillars
    ax.plot3D([-1, -1], [-1, -1], [-1, 1], 'k--')
    ax.plot3D([1, 1], [-1, -1], [-1, 1], 'k--')
    ax.plot3D([-1, -1], [1, 1], [-1, 1], 'k--')
    ax.plot3D([1, 1], [1, 1], [-1, 1], 'k--')

    # guiding lines
    ax.plot3D([-1, 1], [-1, -1], [0, 0], 'k--', alpha=0.4)
    ax.plot3D([-1, 1], [1, 1], [0, 0], 'k--', alpha=0.4)
    ax.plot3D([-1, -1], [-1, 1], [0, 0], 'k--', alpha=0.4)
    ax.plot3D([1, 1], [-1, 1], [0, 0], 'k--', alpha=0.4)

    ax.plot3D([0, 0], [-1, -1], [-1, 1], 'k--', alpha=0.4)
    ax.plot3D([0, 0], [1, 1], [-1, 1], 'k--', alpha=0.4)
    ax.plot3D([0, 0], [-1, 1], [-1, -1], 'k--', alpha=0.4)
    ax.plot3D([0, 0], [1, -1], [1, 1], 'k--', alpha=0.4)

    ax.plot3D([-1, -1], [0, 0], [-1, 1], 'k--', alpha=0.4)
    ax.plot3D([1, 1], [0, 0], [-1, 1], 'k--', alpha=0.4)
    ax.plot3D([-1, 1], [0, 0], [-1, -1], 'k--', alpha=0.4)
    ax.plot3D([-1, 1], [0, 0], [1, 1], 'k--', alpha=0.4)

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
    ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.set_zticks([-1, 0, 1])

    ax.azim = -56
    ax.elev = 15

    origin = (0, 0, 0)
    add_scatter_point(ax, *origin, 'green')

    return ax

def full_factorial():
    designs = doe_fullfact([3] * 3) - 1
    designs = [list(di) for di in designs]
    corners = filter(lambda ci: ci.count(0)==0, designs)
    edges = filter(lambda ci: ci.count(0)==1, designs)
    planecenter = filter(lambda ci: ci.count(0)==2, designs)

    ax = get_unit_cube()
    ax.set_title('Full factorial')

    print("Blue")
    for ci in corners:
        print('\t', ci, hypot(*ci))
        add_scatter_point(ax, *ci, 'lightblue')

    print("Yellow")
    for ci in edges:
        print('\t', ci, hypot(*ci))
        add_scatter_point(ax, *ci, 'yellow')

    print("Purple")
    for ci in planecenter:
        print('\t', ci, hypot(*ci))
        add_scatter_point(ax, *ci, 'purple')

    plt.savefig(DIR_MEDIA / "doe_ff.pdf", bbox_inches="tight")
    plt.show()

def boxbehnken():
    designs = doe_bbdesign(3)
    edges = filter(lambda ci: ci.count(0) == 1, designs)
    ax = get_unit_cube()
    ax.set_title('Box-Behnken')
    for ci in edges:
        add_scatter_point(ax, *ci, 'yellow')

    plt.savefig(DIR_MEDIA / "doe_bb.pdf", bbox_inches="tight")
    plt.show()

def central_composite():
    designs = doe_ccf(3)
    corners = filter(lambda ci: ci.count(0)==0, designs)
    planecenter = filter(lambda ci: ci.count(0)==2, designs)

    ax = get_unit_cube()
    ax.set_title('Central composite')

    for ci in corners:
        add_scatter_point(ax, *ci, 'lightblue')

    for ci in planecenter:
        add_scatter_point(ax, *ci, 'purple')

    plt.savefig(DIR_MEDIA / "doe_ccf.pdf", bbox_inches="tight")
    plt.show()

def plackett_burman():
    designs = doe_pbdesign(3)
    designs = [list(di) for di in designs]
    corners = filter(lambda ci: ci.count(0)==0, designs)
    edges = filter(lambda ci: ci.count(0)==1, designs)
    planecenter = filter(lambda ci: ci.count(0)==2, designs)

    ax = get_unit_cube()
    ax.set_title('Plackett-Burman')

    print("Blue")
    for ci in corners:
        print('\t', ci, hypot(*ci))
        add_scatter_point(ax, *ci, 'lightblue')

    plt.savefig(DIR_MEDIA / "doe_pb.pdf", bbox_inches="tight")

if __name__ == "__main__":
    # full_factorial()
    # boxbehnken()
    # central_composite()
    plackett_burman()


