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
domain_values = {"iron": muopt, "vacuum": mu0}
values_list = [domain_values[mat] for mat in mesh.GetMaterials()]
mu = CoefficientFunction(values_list)

# Boundary condition along GammaH
GammaH = "GammaH1|GammaH2|GammaH3|GammaH4"
Hy = 10000
fi0 = -Hy * y

# Init R, H, B inside iron
fes1 = H1(mesh, order=3, definedon="iron")
fes1 = Compress(fes1)
Rx = GridFunction(fes1)
Ry = GridFunction(fes1)
Rx.Set(0.0)
Ry.Set(0.0)
Rxlast = GridFunction(fes1)
Rylast = GridFunction(fes1)
Rxlast.Set(0.0)
Rylast.Set(0.0)
Hx = GridFunction(fes1)
Hy = GridFunction(fes1)
Hx.Set(0.0)
Hy.Set(0.0)
Bx = GridFunction(fes1)
By = GridFunction(fes1)
Bx.Set(0.0)
By.Set(0.0)

fes = H1(mesh, order=3, dirichlet=GammaH)
u = fes.TrialFunction()
v = fes.TestFunction()
gfu = GridFunction(fes)
gfu.Set(fi0, BND)

# A
a = BilinearForm(fes)
a += mu * grad(u) * grad(v) * dx
a.Assemble()

f = LinearForm(fes)
f += (Rx * grad(v)[0] + Ry * grad(v)[1]) * dx("iron")
f.Assemble()

r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

# Nonlinear fixed point iteration
err = 1
rxlast = 0
rylast = 0
while err > 1e-8:
    Hx.Set(-grad(gfu)[0])
    Hy.Set(-grad(gfu)[1])

    Bx.Set(2.0 * Bs / pi * atan(Hx / H0))
    By.Set(2.0 * Bs / pi * atan(Hy / H0))

    Rxlast.Set(Rx)
    Rylast.Set(Ry)

    Rx.Set(Bx - muopt * Hx)
    Ry.Set(By - muopt * Hy)

    err = 0
    for i in range(len(Rx.vec)):
        err += sqrt((Rx.vec.data[i] - Rxlast.vec.data[i]) ** 2 + (Ry.vec.data[i] - Rylast.vec.data[i]) ** 2)
    err /= len(Rx.vec)
    print(err)

    f = LinearForm(fes)
    f += (Rx * grad(v)[0] + Ry * grad(v)[1]) * dx("iron")
    f.Assemble()

    r = f.vec.CreateVector()
    r.data = f.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

Draw((Hx, Hy), mesh, "H")
Draw((Bx, By), mesh, "B")
Draw(gfu, mesh, "H1")
