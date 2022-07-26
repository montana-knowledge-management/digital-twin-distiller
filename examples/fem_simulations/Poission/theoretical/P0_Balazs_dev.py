# PREAMBLE

from netgen.geom2d import SplineGeometry
from ngsolve import *

geo = SplineGeometry()


# METADATA

# empty


# MATERIAL DEFINITIONS

# empty


# BOUNDARY DEFINITIONS

# empty


# GEOMETRY

geo.AppendPoint(0, 0)
geo.AppendPoint(1, 0)
geo.AppendPoint(1, 1)
geo.AppendPoint(0, 1)
geo.Append(["line", 0, 1], leftdomain=1, rightdomain=0, bc="1")
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0, bc="2")
geo.Append(["line", 2, 3], leftdomain=1, rightdomain=0, bc="3")
geo.Append(["line", 3, 0], leftdomain=1, rightdomain=0, bc="4")
mesh = Mesh(geo.GenerateMesh(maxh=0.025))

fes = H1(mesh, order=1, dirichlet=[1, 2, 3, 4])
gfu = GridFunction(fes)

alpha = 60
r0 = 1
x0 = 1.25
y0 = -0.25
g = atan(alpha * (sqrt((x - x0) ** 2 + (y - y0) ** 2) - r0))
gfu.Set(g, BND)

u = fes.TrialFunction()
v = fes.TestFunction()
a = BilinearForm(fes, symmetric=True)
a += grad(u) * grad(v) * dx
a.Assemble()

f = LinearForm(fes)
f += (
    -(
        alpha * (alpha**2 * (sqrt((x - x0) ** 2 + (y - y0) ** 2) - r0) ** 2 + 1)
        - 2 * alpha**3 * sqrt((x - x0) ** 2 + (y - y0) ** 2) * (sqrt((x - x0) ** 2 + (y - y0) ** 2) - r0)
    )
    / (sqrt((x - x0) ** 2 + (y - y0) ** 2) * (alpha**2 * (sqrt((x - x0) ** 2 + (y - y0) ** 2) - r0) ** 2 + 1) ** 2)
    * v
    * dx
)
# f +=  ((alpha*(y-y0))/(sqrt((y-y0)**2+(x-x0)**2)*(alpha**2*(sqrt((y-y0)**2+(x-x0)**2)-r0)**2+1)))*v*ds(definedon="3")
# f += -((alpha*(x-x0))/(sqrt((x-x0)**2+(y-y0)**2)*(alpha**2*(sqrt((x-x0)**2+(y-y0)**2)-r0)**2+1)))*v*ds(definedon="4")
f.Assemble()

r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r
Draw(gfu)
Draw(grad(gfu), mesh, "grad")
# BLOCK LABELS

# empty


# SOLVER

# fes = H1(mesh, order=3, dirichlet="gnd")
# u = fes.TrialFunction()
# v = fes.TestFunction()
# f = LinearForm(fes)
# f += 32 * (y*(1-y)+x*(1-x)) * v * dx
# a = BilinearForm(fes, symmetric=True)
# a += grad(u)*grad(v)*dx
# a.Assemble()
# f.Assemble()
# gfu = GridFunction(fes)
# gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
# Draw (gfu)


# POSTPROCESSING AND EXPORTING

# empty


# CLOSING STEPS

# empty
