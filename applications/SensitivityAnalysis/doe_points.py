import matplotlib.pyplot as plt
import itertools
from model import *
from math import hypot

def add_scatter_point(ax, x, y, z, color):
    ax.scatter(x, y, z, c=color, s=70, edgecolor='k', alpha=1)

def get_unit_cube():
    fig = plt.figure()
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
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    ax.set_zticks([-1, 0, 1])

    ax.azim = -56
    ax.elev = 15
    return ax

def full_factorial():
    N = 3
    parameters = [[-1, 0, 1]]* N
    cases = list(itertools.product(*parameters))

    # partitioning
    origin = (0, 0, 0)
    corners = filter(lambda ci: ci.count(0)==0, cases)
    edges = filter(lambda ci: ci.count(0)==1, cases)
    planecenter = filter(lambda ci: ci.count(0)==2, cases)

    ax = get_unit_cube()
    ax.set_title('Full factorial')
    add_scatter_point(ax, *origin, 'green')

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

    plt.savefig(DIR_MEDIA / "doe_full.png", bbox_inches="tight")
    plt.show()

def boxbehnken():
    N = 3
    parameters = [[-1, 0, 1]]* N
    cases = list(itertools.product(*parameters))

    # partitioning
    origin = (0, 0, 0)
    edges = filter(lambda ci: ci.count(0)==1, cases)

    ax = get_unit_cube()
    ax.set_title('Box Behnken')
    add_scatter_point(ax, *origin, 'green')

    for ci in edges:
        add_scatter_point(ax, *ci, 'yellow')


    plt.savefig(DIR_MEDIA / "doe_boxbehnken.png", bbox_inches="tight")
    plt.show()

def koshal():
    N = 3
    parameters = [[-1, 0, 1]]* N
    cases = list(itertools.product(*parameters))

    # partitioning
    origin = (0, 0, 0)
    edges = filter(lambda ci: ci.count(0)==1, cases)
    edges = filter(lambda ci: ci.count(1)==N-1, edges)
    planecenter = filter(lambda ci: ci.count(0)==2, cases)


    ax = get_unit_cube()
    ax.set_title('Koshal')
    add_scatter_point(ax, *origin, 'green')

    for ci in edges:
        add_scatter_point(ax, *ci, 'yellow')

    for ci in planecenter:
        add_scatter_point(ax, *ci, 'purple')

    plt.savefig(DIR_MEDIA / "doe_koshal.png", bbox_inches="tight")
    plt.show()

def central_composite():
    N = 3
    parameters = [[-1, 0, 1]]* N
    cases = list(itertools.product(*parameters))

    # partitioning
    origin = (0, 0, 0)
    corners = filter(lambda ci: ci.count(0)==0, cases)
    planecenter = filter(lambda ci: ci.count(0)==2, cases)

    ax = get_unit_cube()
    ax.set_title('Central composite')
    add_scatter_point(ax, *origin, 'green')

    for ci in corners:
        add_scatter_point(ax, *ci, 'lightblue')

    for ci in planecenter:
        add_scatter_point(ax, *ci, 'purple')

    plt.savefig(DIR_MEDIA / "doe_centralcomposite.png", bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    full_factorial()
    # boxbehnken()
    # koshal()
    # central_composite()

