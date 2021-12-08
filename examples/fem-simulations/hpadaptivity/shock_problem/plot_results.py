from digital_twin_distiller import ModelDir
from numpy.core.function_base import linspace
from digital_twin_distiller import setup_matplotlib
import matplotlib.pyplot as plt
import numpy as np

ModelDir.set_base(__file__)

setup_matplotlib()

x = linspace(0, 1, 1001)
y = linspace(0, 1, 1001)

x_0 = np.zeros_like(x)
x_1 = np.ones_like(x)

y_0 = np.zeros_like(y)
y_1 = np.ones_like(y)

xx, yy = np.meshgrid(x, y)
alpha = 60
r0 = np.array([1.25, -0.25])

u = lambda xi, yi: np.arctan(alpha * (np.sqrt((r0[0]-xi)**2 + (r0[1]-yi)**2)-1))

# plt.figure(figsize=(10, 10))
# plt.contourf(xx, yy, u(xx, yy))
# plt.xlim(0, 1)
# plt.ylim(0, 1)


plt.show()
plt.figure(figsize=(9, 9))

plt.subplot(335)
plt.contourf(xx, yy, u(xx, yy))
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.axis('off')

plt.subplot(338)
plt.plot(x, u(x, y_0))
plt.xlabel('x')
plt.ylabel(r'$\Gamma_1$')

plt.subplot(332)
plt.plot(x, u(x, y_1))
plt.xlabel('x')
plt.ylabel(r'$\Gamma_3$')

plt.subplot(334)
plt.plot(-u(x_0, y), y)
plt.ylabel('y')
plt.xlabel(r'$\Gamma_4$')

plt.subplot(336)
plt.plot(u(x_1, y), y)
plt.ylabel('y')
plt.xlabel(r'$\Gamma_2$')

plt.tight_layout()
plt.savefig(ModelDir.MEDIA / "boundaries.png", dpi=450, bbox_inches="tight")
plt.show()
