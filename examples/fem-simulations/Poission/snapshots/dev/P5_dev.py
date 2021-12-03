from ngsolve import *
from netgen.geom2d import SplineGeometry

# Geometry
geo = SplineGeometry()
geo.AppendPoint(-2, -2)
geo.AppendPoint(2, -2)
geo.AppendPoint(2,  2)
geo.AppendPoint(-2,  2)
geo.AppendPoint(-0.5,  -0.5)
geo.AppendPoint(0.5,  -0.5)
geo.AppendPoint(0.5, 0.5)
geo.AppendPoint(-0.5, 0.5)
geo.Append(["line", 0, 1], leftdomain=1, rightdomain=0)
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0)
geo.Append(["line", 2, 3], leftdomain=1, rightdomain=0)
geo.Append(["line", 3, 0], leftdomain=1, rightdomain=0)
geo.Append(["line", 4, 5], leftdomain=2, rightdomain=1)
geo.Append(["line", 5, 6], leftdomain=2, rightdomain=1)
geo.Append(["line", 6, 7], leftdomain=2, rightdomain=1)
geo.Append(["line", 7, 4], leftdomain=2, rightdomain=1)
geo.SetMaterial(1, "outer")
geo.SetMaterial(2, "inner")
mesh = Mesh(geo.GenerateMesh(maxh=0.05))

# Constants
domain_values = {'inner': 1000*4*pi*1e-7,  'outer': 4*pi*1e-7}
values_list   = [domain_values[mat] for mat in mesh.GetMaterials()]
mu            = CoefficientFunction(values_list)
H0            = CoefficientFunction([0,-100,0,100,0,0,0,0])
# Mesh
fes = H1(mesh, order=3, dirichlet=[])
u   = fes.TrialFunction()
v   = fes.TestFunction()
gfu = GridFunction(fes)

# A
a = BilinearForm(fes)
a += 1/mu*grad(u)*grad(v)*dx
a.Assemble()

# f
f = LinearForm(fes)
f += H0 * v * ds
f.Assemble()

gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * f.vec
Draw(gfu)