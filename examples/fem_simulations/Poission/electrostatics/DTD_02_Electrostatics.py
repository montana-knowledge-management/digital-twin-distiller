from ngsolve import *
from netgen.geom2d import SplineGeometry

# Geometry
geo = SplineGeometry()
geo.AppendPoint(0, 0)
geo.AppendPoint(0.5, 0)
geo.AppendPoint(2,  0)
geo.AppendPoint(2,  2)
geo.AppendPoint(0,  2)
geo.AppendPoint(0,  0.5)
geo.AppendPoint(0.5, 0.5)
geo.Append(["line", 0, 1], leftdomain=2, rightdomain=0, bc="GammaD1")
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0, bc="GammaD2")
geo.Append(["line", 2, 3], leftdomain=1, rightdomain=0, bc="GammaN1")
geo.Append(["line", 3, 4], leftdomain=1, rightdomain=0, bc="GammaD3")
geo.Append(["line", 4, 5], leftdomain=1, rightdomain=0, bc="GammaN2")
geo.Append(["line", 5, 0], leftdomain=2, rightdomain=0, bc="GammaN3")
geo.Append(["line", 5, 6], leftdomain=1, rightdomain=2)
geo.Append(["line", 6, 1], leftdomain=1, rightdomain=2)
geo.SetMaterial(1, "outer")
geo.SetMaterial(2, "inner")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

# Constant epsilon
eps0          = 8.854e-12
epsr          = 10
domain_values = {'inner': eps0*epsr,  'outer': eps0}
values_list   = [domain_values[mat] for mat in mesh.GetMaterials()]
epsilon       = CoefficientFunction(values_list)

# Constant rho
domain_values = {'inner': 0,  'outer': 0}
values_list   = [domain_values[mat] for mat in mesh.GetMaterials()]
rho           = CoefficientFunction(values_list)

# Dirichlet boundary condition
GammaD = "GammaD1|GammaD2|GammaD3"
g      = 100*(y+2)/4

# Neumann boundary condition
h1  = 0*eps0
h2  = 0*eps0
h3  = 0*eps0

# Mesh
fes = H1(mesh, order=3, dirichlet=GammaD)
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
f += v*rho*dx
f += v*h1*ds(definedon="GammaN1")
f += v*h2*ds(definedon="GammaN2")
f += v*h3*ds(definedon="GammaN3")
f.Assemble()

# Dirichlet set & solve
r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r

# Draw u & E
Draw(gfu)
Draw(-grad(gfu),mesh,"grad")