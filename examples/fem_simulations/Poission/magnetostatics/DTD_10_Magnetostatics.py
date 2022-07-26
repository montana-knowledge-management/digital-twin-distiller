from netgen.geom2d import SplineGeometry
from ngsolve import *

# Geometry
geo = SplineGeometry()
geo.AppendPoint(-2, -2)
geo.AppendPoint(2, -2)
geo.AppendPoint(2, 2)
geo.AppendPoint(-2, 2)
geo.AppendPoint(-0.5, -0.5)
geo.AppendPoint(0.5, -0.5)
geo.AppendPoint(0.5, 0.5)
geo.AppendPoint(-0.5, 0.5)
geo.Append(["line", 0, 1], leftdomain=1, rightdomain=0, bc="GammaH1")
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0, bc="GammaH2")
geo.Append(["line", 2, 3], leftdomain=1, rightdomain=0, bc="GammaH3")
geo.Append(["line", 3, 0], leftdomain=1, rightdomain=0, bc="GammaH4")
geo.Append(["line", 4, 5], leftdomain=2, rightdomain=1)
geo.Append(["line", 5, 6], leftdomain=2, rightdomain=1)
geo.Append(["line", 6, 7], leftdomain=2, rightdomain=1)
geo.Append(["line", 7, 4], leftdomain=2, rightdomain=1)
geo.SetMaterial(1, "vacuum")
geo.SetMaterial(2, "iron")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

# Constant mu
mu0 = 4 * pi * 1e-7
mur = 1000
domain_values = {"iron": mu0 * mur, "vacuum": mu0}
values_list = [domain_values[mat] for mat in mesh.GetMaterials()]
mu = CoefficientFunction(values_list)

# Boundary condition along GammaH
Hx = 100
Hy = 200
K = CoefficientFunction([-Hx, -Hy, Hx, Hy])

# Mesh
fes = H1(mesh, order=3, dirichlet=[])
u = fes.TrialFunction()
v = fes.TestFunction()
gfu = GridFunction(fes)

# A
a = BilinearForm(fes)
a += 1 / mu * grad(u) * grad(v) * dx
a.Assemble()

# f
f = LinearForm(fes)
f += v * K * ds
f.Assemble()

gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
Draw(gfu)
Draw((grad(gfu)[1], -grad(gfu)[0]), mesh, "B")
