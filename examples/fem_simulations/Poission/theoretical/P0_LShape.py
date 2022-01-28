# PREAMBLE

from ngsolve import *
from netgen.geom2d import SplineGeometry

dudx = 2/3*(x*sin(2/3*atan2(y,x)+pi/3)-y*cos(2/3*atan2(y,x)+pi/3))/((x**2+y**2)**(2/3))
dudy = 2/3*(y*sin(2/3*atan2(y,x)+pi/3)+x*cos(2/3*atan2(y,x)+pi/3))/((x**2+y**2)**(2/3))

geo = SplineGeometry()
geo.AppendPoint( 0,  0)
geo.AppendPoint( 0, -1)
geo.AppendPoint( 1, -1)
geo.AppendPoint( 1,  1)
geo.AppendPoint(-1,  1)
geo.AppendPoint(-1,  0)

geo.Append(["line", 0, 1], leftdomain=1, rightdomain=0, bc="GammaD1")
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0, bc="GammaN1")
geo.Append(["line", 2, 3], leftdomain=1, rightdomain=0, bc="GammaN2")
geo.Append(["line", 3, 4], leftdomain=1, rightdomain=0, bc="GammaN3")
geo.Append(["line", 4, 5], leftdomain=1, rightdomain=0, bc="GammaN4")
geo.Append(["line", 5, 0], leftdomain=1, rightdomain=0, bc="GammaD2")
mesh = Mesh(geo.GenerateMesh(maxh=0.025))

fes = H1(mesh, order=10, dirichlet="GammaD1|GammaD2")
gfu = GridFunction(fes)
gfu.Set(0, BND)

u = fes.TrialFunction()
v = fes.TestFunction()
a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx
a.Assemble()

f = LinearForm(fes)
f += -dudy*v*ds(definedon="GammaN1")
f +=  dudy*v*ds(definedon="GammaN3")
f +=  dudx*v*ds(definedon="GammaN2")
f += -dudx*v*ds(definedon="GammaN4")
f.Assemble()

r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r
Draw(gfu)


print(mesh.Elements(VOL))
print(len(gfu.vec.data))
print(gfu(1e-3,1e-3))
print(a.Energy(gfu.vec))