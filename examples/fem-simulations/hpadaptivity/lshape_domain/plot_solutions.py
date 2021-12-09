from digital_twin_distiller import setup_matplotlib
import matplotlib.pyplot as plt
import numpy as np


x = np.linspace(-1, 1, 1000)
y = np.linspace(-1, 1, 1000)

xx, yy = np.meshgrid(x, y)
# idx3 = np.logical_or(xx <= 0.0, yy <=0.0)
# idx3 = np.logical_not(idx3)

# xx = xx[idx3]
# yy = yy[idx3]

plt.figure(figsize=(10, 10))

r = np.hypot(xx, yy)
phi = np.arctan2(yy, xx)
u = r**(2/3) * np.sin(2/3*phi+np.pi/3)
plt.contourf(xx, yy, u)

plt.colorbar()
plt.show()


