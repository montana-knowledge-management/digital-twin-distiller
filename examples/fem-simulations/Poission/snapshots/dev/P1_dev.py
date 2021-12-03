# PREAMBLE

from ngsolve import *
from netgen.geom2d import SplineGeometry

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
geo.Append(["line", 0, 1], leftdomain=1, rightdomain=0)
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0)
geo.Append(["line", 2, 3], leftdomain=1, rightdomain=0)
geo.Append(["line", 3, 0], leftdomain=1, rightdomain=0)
mesh = Mesh(geo.GenerateMesh(maxh=0.05))

fes = H1(mesh, order=2, dirichlet=[1,3,4])
gfu = GridFunction(fes)

g = y*(1-y)
gfu.Set(g, BND)

u = fes.TrialFunction()
v = fes.TestFunction()
a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx
a.Assemble()
f = LinearForm(fes)
f += 1*v*dx + 1*v*ds
f.Assemble()

r = f.vec.CreateVector()
r.data = f.vec - a.mat * gfu.vec

gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * r
Draw(gfu)

# BLOCK LABELS

# empty



# SOLVER

#fes = H1(mesh, order=3, dirichlet="gnd")
#u = fes.TrialFunction()
#v = fes.TestFunction()
#f = LinearForm(fes)
#f += 32 * (y*(1-y)+x*(1-x)) * v * dx
#a = BilinearForm(fes, symmetric=True)
#a += grad(u)*grad(v)*dx
#a.Assemble()
#f.Assemble()
#gfu = GridFunction(fes)
#gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
#Draw (gfu)


# POSTPROCESSING AND EXPORTING

# empty



# CLOSING STEPS

# empty



