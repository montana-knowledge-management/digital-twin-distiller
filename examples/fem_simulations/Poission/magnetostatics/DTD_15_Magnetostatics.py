from netgen.geom2d import SplineGeometry
from ngsolve import *

# Geometry
geo = SplineGeometry()
geo.AppendPoint(0, 0)
geo.AppendPoint(0.5, 0)
geo.AppendPoint(2, 0)
geo.AppendPoint(2, 2)
geo.AppendPoint(0, 2)
geo.AppendPoint(0, 0.5)
geo.AppendPoint(0.5, 0.5)
geo.Append(["line", 0, 1], leftdomain=2, rightdomain=0, bc="GammaH1")
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0, bc="GammaH2")
geo.Append(["line", 2, 3], leftdomain=1, rightdomain=0, bc="GammaH3")
geo.Append(["line", 3, 4], leftdomain=1, rightdomain=0, bc="GammaH4")
geo.Append(["line", 4, 5], leftdomain=1, rightdomain=0, bc="GammaB1")
geo.Append(["line", 5, 0], leftdomain=2, rightdomain=0, bc="GammaB2")
geo.Append(["line", 5, 6], leftdomain=1, rightdomain=2)
geo.Append(["line", 6, 1], leftdomain=1, rightdomain=2)
geo.SetMaterial(1, "vacuum")
geo.SetMaterial(2, "iron")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

# Constants
mu0 = 4 * pi * 1e-7
Bs = 2
H0 = 100
mumax = 2 * Bs / pi / H0
muopt = 1.1 * mumax / 2
nuopt = 1 / muopt
domain_values = {"iron": nuopt, "vacuum": 1 / mu0}
values_list = [domain_values[mat] for mat in mesh.GetMaterials()]
nu = CoefficientFunction(values_list)

# Boundary condition
GammaB = "GammaB1|GammaB2"
Hy = 100
K = CoefficientFunction([0, 0, -Hy, 0])

# Init I, H, B inside iron
fes1 = H1(mesh, order=3, definedon="iron")
fes1 = Compress(fes1)
Ix = GridFunction(fes1)
Iy = GridFunction(fes1)
Ix.Set(0.0)
Iy.Set(0.0)
Ixlast = GridFunction(fes1)
Iylast = GridFunction(fes1)
Ixlast.Set(0.0)
Iylast.Set(0.0)
Hx = GridFunction(fes1)
Hy = GridFunction(fes1)
Hx.Set(0.0)
Hy.Set(0.0)
Bx = GridFunction(fes1)
By = GridFunction(fes1)
Bx.Set(0.0)
By.Set(0.0)

fes = H1(mesh, order=3, dirichlet=GammaB)
u = fes.TrialFunction()
v = fes.TestFunction()
gfu = GridFunction(fes)
gfu.Set(0, BND)

# A
a = BilinearForm(fes)
a += nu * grad(u) * grad(v) * dx
a.Assemble()

f = LinearForm(fes)
f += K * v * ds
f += -(Ix * grad(v)[1] - Iy * grad(v)[0]) * dx("iron")
f.Assemble()

r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

# Nonlinear fixed point iteration
err = 1
ixlast = 0
iylast = 0
while err > 1e-8:
    Bx.Set(grad(gfu)[1])
    By.Set(-grad(gfu)[0])

    Hx.Set(nuopt * Bx + Ix)
    Hy.Set(nuopt * By + Iy)

    Ixlast.Set(Ix)
    Iylast.Set(Iy)

    Ix.Set(Hx - nuopt * (2.0 * Bs / pi * atan(Hx / H0)))
    Iy.Set(Hy - nuopt * (2.0 * Bs / pi * atan(Hy / H0)))

    err = 0
    for i in range(len(Ix.vec)):
        err += sqrt((Ix.vec.data[i] - Ixlast.vec.data[i]) ** 2 + (Iy.vec.data[i] - Iylast.vec.data[i]) ** 2)
    err /= len(Ix.vec)
    print(err)

    f = LinearForm(fes)
    f += K * v * ds
    f += -(Ix * grad(v)[1] - Iy * grad(v)[0]) * dx("iron")
    f.Assemble()

    r = f.vec.CreateVector()
    r.data = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

Draw((Hx, Hy), mesh, "H")
Draw((Bx, By), mesh, "B")
Draw((grad(gfu)[1], -grad(gfu)[0]), mesh, "B1")
Draw(sqrt((grad(gfu)[1]) ** 2 + (grad(gfu)[0]) ** 2), mesh, "abs(B)")
