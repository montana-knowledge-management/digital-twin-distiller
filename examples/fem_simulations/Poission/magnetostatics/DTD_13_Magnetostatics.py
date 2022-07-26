from netgen.geom2d import SplineGeometry
from ngsolve import *

# Geometry
geo = SplineGeometry()
geo.AppendPoint(0, 0)
geo.AppendPoint(0.5, 0)
geo.AppendPoint(1, 0)
geo.AppendPoint(1.1, 0)
geo.AppendPoint(2, 0)
geo.AppendPoint(2, 2)
geo.AppendPoint(0, 2)
geo.AppendPoint(0, 0.5)
geo.AppendPoint(0.5, 0.5)
geo.AppendPoint(1, 0.05)
geo.AppendPoint(1.1, 0.05)
geo.Append(["line", 0, 1], leftdomain=2, rightdomain=0, bc="GammaH1")
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0, bc="GammaH2")
geo.Append(["line", 2, 3], leftdomain=3, rightdomain=0, bc="GammaH3")
geo.Append(["line", 3, 4], leftdomain=1, rightdomain=0, bc="GammaH4")
geo.Append(["line", 4, 5], leftdomain=1, rightdomain=0, bc="GammaB1")
geo.Append(["line", 5, 6], leftdomain=1, rightdomain=0, bc="GammaB2")
geo.Append(["line", 6, 7], leftdomain=1, rightdomain=0, bc="GammaB3")
geo.Append(["line", 7, 0], leftdomain=2, rightdomain=0, bc="GammaB4")
geo.Append(["line", 7, 8], leftdomain=1, rightdomain=2)
geo.Append(["line", 8, 1], leftdomain=1, rightdomain=2)
geo.Append(["line", 2, 9], leftdomain=1, rightdomain=3)
geo.Append(["line", 9, 10], leftdomain=1, rightdomain=3)
geo.Append(["line", 3, 10], leftdomain=3, rightdomain=1)
geo.SetMaterial(1, "air")
geo.SetMaterial(2, "iron")
geo.SetMaterial(3, "coil")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

# Constant mu
mu0 = 4 * pi * 1e-7
mur = 1000
domain_values = {"iron": mu0 * mur, "air": mu0, "coil": mu0}
values_list = [domain_values[mat] for mat in mesh.GetMaterials()]
mu = CoefficientFunction(values_list)

# Constant J0
domain_values = {"iron": 0, "air": 0, "coil": 1e6}
values_list = [domain_values[mat] for mat in mesh.GetMaterials()]
J0 = CoefficientFunction(values_list)

# Boundary conditions
GammaB = "GammaB1|gammaB2|GammaB3|GammaB4"
K = CoefficientFunction([0, 0, 0, 0])

# Mesh
fes = H1(mesh, order=2, dirichlet=GammaB)
u = fes.TrialFunction()
v = fes.TestFunction()
gfu = GridFunction(fes)
gfu.Set(0, BND)

R0 = GridFunction(fes)
R0.Set(x)

# A
a = BilinearForm(fes)
a += 1 / R0 * 1 / mu * grad(u) * grad(v) * dx
a.Assemble()

# f
f = LinearForm(fes)
f += v * J0 * dx + v * K * ds
f.Assemble()

# Dirichlet set & solve
r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r
for i in range(len(gfu.vec.data)):
    x = mesh.vertices[i].point[0]
    if x > 1e-6:
        gfu.vec.data[i] = gfu.vec.data[i] / x

Draw(gfu)
Draw((grad(gfu)[1], -grad(gfu)[0]), mesh, "B")
