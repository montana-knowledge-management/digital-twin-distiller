from ngsolve import *
from netgen.geom2d import SplineGeometry

# Geometry
geo = SplineGeometry()
geo.AppendPoint(0, -2)
geo.AppendPoint(2, -2)
geo.AppendPoint(2,  2)
geo.AppendPoint(0,  2)
geo.AppendPoint(0,  0.5)
geo.AppendPoint(0.5,  0.5)
geo.AppendPoint(0.5, -0.5)
geo.AppendPoint(0, -0.5)
geo.Append(["line", 0, 1], leftdomain=1, rightdomain=0)
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0)
geo.Append(["line", 2, 3], leftdomain=1, rightdomain=0)
geo.Append(["line", 3, 4], leftdomain=1, rightdomain=0)
geo.Append(["line", 4, 7], leftdomain=2, rightdomain=0)
geo.Append(["line", 7, 0], leftdomain=1, rightdomain=0)
geo.Append(["line", 4, 5], leftdomain=1, rightdomain=2)
geo.Append(["line", 5, 6], leftdomain=1, rightdomain=2)
geo.Append(["line", 6, 7], leftdomain=1, rightdomain=2)
geo.SetMaterial(1, "outer")
geo.SetMaterial(2, "inner")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

# Constants
domain_values = {'inner': 8.854e-12*10,  'outer': 8.854e-12}
values_list   = [domain_values[mat] for mat in mesh.GetMaterials()]
epsilon       = CoefficientFunction(values_list)

g   = 100*(y+2)/4 # Dirichlet

rho = 0

# Mesh
fes = H1(mesh, order=3, dirichlet=[1,3])
u   = fes.TrialFunction()
v   = fes.TestFunction()
gfu = GridFunction(fes)
gfu.Set(g, BND)

# A
a = BilinearForm(fes)
a += epsilon*grad(u)*grad(v)*dx
a.Assemble()

# f
f = LinearForm(fes)
f += rho * v * dx
f.Assemble()

# Dirichlet set & solve
r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r
Draw(gfu)