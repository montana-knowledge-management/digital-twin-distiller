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

# Constant mu
mu0           = 4*pi*1e-7
mur           = 1000
domain_values = {'iron': mu0*mur, 'vacuum': mu0}
values_list   = [domain_values[mat] for mat in mesh.GetMaterials()]
mu            = CoefficientFunction(values_list)

# Boundary condition along GammaH
Hy = 200
K  = CoefficientFunction([0,0,-Hy,0])

GammaB = "GammaB1|GammaB2"

# Mesh
fes = H1(mesh, order=3, dirichlet=GammaB)
u   = fes.TrialFunction()
v   = fes.TestFunction()
gfu = GridFunction(fes)
gfu.Set(0, BND)

# A
a = BilinearForm(fes)
a += 1/mu*grad(u)*grad(v)*dx
a.Assemble()

# f
f = LinearForm(fes)
f += v*K*ds
f.Assemble()

r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r
Draw(gfu)
Draw((grad(gfu)[1],-grad(gfu)[0]),mesh,"B")