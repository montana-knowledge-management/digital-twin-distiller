import numpy as np
import matplotlib.pyplot as plt
from digital_twin_distiller import ModelDir
from scipy import interpolate


ModelDir.set_base(__file__)

x0 = 0.5
y0 = 0.5
alpha = 50
r0 = 0.25

def u(x, y):
    return np.arctan(alpha*(np.hypot(x-x0, y-y0)-r0))


def get_core_points(filename, nb_x=250, nb_y=250):
    x, y, V = np.genfromtxt(filename, unpack=True, delimiter=',')
    xi = np.linspace(x.min(), x.max(), nb_x)
    yi = np.linspace(y.min(), y.max(), nb_y)
    zi = interpolate.griddata((x, y), V, (xi[None, :], yi[:, None]), method='cubic')

    xx, yy = np.meshgrid(xi, yi)

    return xx, yy, zi

xx, yy, zz_agros = get_core_points(ModelDir.DATA/'r_agros.csv')
xx, yy, zz_femm = get_core_points(ModelDir.DATA/'r_femm.csv')

plt.figure(figsize=(10, 5))

plt.subplot(231)
plt.contourf(xx, yy, zz_agros, cmap='inferno')
plt.colorbar()
plt.axis('off')
plt.title("Agros2D")

plt.subplot(232)
plt.contourf(xx, yy, zz_femm, cmap='inferno')
plt.colorbar()
plt.axis('off')
plt.title("FEMM")

plt.subplot(233)
plt.contourf(xx, yy, u(xx, yy), cmap='inferno')
plt.colorbar()
plt.axis('off')
plt.title("Exact")

plt.subplot(234)
plt.contourf(xx, yy, u(xx, yy)-zz_agros, cmap='inferno')
plt.colorbar()
plt.axis('off')
plt.title("Difference: uexact-agros")

plt.subplot(235)
plt.contourf(xx, yy, u(xx, yy)-zz_femm, cmap='inferno')
plt.colorbar()
plt.axis('off')
plt.title("Difference: uexact-femm")

plt.subplot(236)
plt.contourf(xx, yy, zz_agros-zz_femm, cmap='inferno')
plt.colorbar()
plt.axis('off')
plt.title("Difference: agros-femm")

plt.tight_layout()
plt.savefig(ModelDir.MEDIA / "diff.png", dpi=450, bbox_inches="tight")
plt.show()


"""
exact,
femm,
agros,
exact-agros,
exact-femm,
agros-femm
"""


"""
G1: y=0.5 atan(50*(sqrt((x-0.5)**2)-0.25))
G2: x=1 atan(50*(sqrt(0.5**2+(y-0.5)**2)-0.25))
G3: y=1 atan(50*(sqrt((x-0.5)**2+0.5**2)-0.25))
G4: x=0.5 atan(50*(sqrt((y-0.5)**2)-0.25))
"""
