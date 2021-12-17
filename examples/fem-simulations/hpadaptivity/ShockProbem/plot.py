import numpy as np
import matplotlib.pyplot as plt
from math import atan, hypot, sqrt



x = np.linspace(0.5, 1, 1001)
y = np.linspace(0.5, 1, 1001)

x0 = 0.5
y0 = 0.5
alpha = 50
r0 = 0.25

def u(x, y):
    return np.arctan(alpha*(np.hypot(x-x0, y-y0)-r0))


def u1(x, y):
    return atan(alpha*(hypot(x-x0, y-y0)-r0))

def u3(x,y):
    return atan(alpha*(sqrt((x-x0)**2+(y-y0)**2)-r0))

# xx, yy = np.meshgrid(x, y)

# u = np.arctan(alpha*(np.hypot(xx-x0, yy-y0)-r0))

# plt.figure()
# plt.contourf(xx, yy, u, cmap='jet')
# plt.colorbar()
# plt.tight_layout()
# plt.show()
u4 = lambda x,y: atan(alpha*(sqrt((y-y0)**2)-r0))

u_ = [u(xi, 0.5) for xi in x ]
u_3 = [u3(xi, 0.5) for xi in x]

plt.figure()
plt.plot(y, u_3)
plt.show()


"""
G1: y=0.5 atan(50*(sqrt((x-0.5)**2)-0.25))
G2: x=1 atan(50*(sqrt(0.5**2+(y-0.5)**2)-0.25))
G3: y=1 atan(50*(sqrt((x-0.5)**2+0.5**2)-0.25))
G4: x=0.5 atan(50*(sqrt((y-0.5)**2)-0.25))
"""
