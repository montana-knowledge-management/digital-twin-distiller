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

geo.AppendPoint(8.79, 4.35)
geo.AppendPoint(10.39, 6.69)
geo.AppendPoint(6.38, 5.57)
geo.AppendPoint(11.33, 0.23)
geo.AppendPoint(13.75, 2.34)
geo.AppendPoint(12.62, 3.34)
geo.AppendPoint(6.57, 1.22)
geo.AppendPoint(8.04, 3.3)
geo.AppendPoint(6.41, 4.16)
geo.AppendPoint(4.3, 3.68)
geo.AppendPoint(14.03, 6.09)
geo.Append(["line", 0, 1], leftdomain=1, rightdomain=3, bc="gnd")
geo.Append(["line", 0, 7], leftdomain=3, rightdomain=2, bc="gnd")
geo.Append(["line", 1, 2], leftdomain=1, rightdomain=0, bc="gnd")
geo.Append(["line", 1, 10], leftdomain=0, rightdomain=5, bc="gnd")
geo.Append(["line", 2, 0], leftdomain=1, rightdomain=2, bc="gnd")
geo.Append(["line", 2, 9], leftdomain=2, rightdomain=0, bc="gnd")
geo.Append(["line", 3, 4], leftdomain=4, rightdomain=0, bc="gnd")
geo.Append(["line", 4, 5], leftdomain=4, rightdomain=5, bc="gnd")
geo.Append(["line", 5, 6], leftdomain=4, rightdomain=3, bc="gnd")
geo.Append(["line", 5, 1], leftdomain=3, rightdomain=5, bc="gnd")
geo.Append(["line", 6, 3], leftdomain=4, rightdomain=0, bc="gnd")
geo.Append(["line", 7, 8], leftdomain=3, rightdomain=2, bc="gnd")
geo.Append(["line", 8, 9], leftdomain=3, rightdomain=2, bc="gnd")
geo.Append(["line", 9, 6], leftdomain=3, rightdomain=0, bc="gnd")
geo.Append(["line", 10, 4], leftdomain=0, rightdomain=5, bc="gnd")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))


# BLOCK LABELS

# empty



# SOLVER

fes = H1(mesh, order=3, dirichlet="gnd")
u = fes.TrialFunction()
v = fes.TestFunction()
f = LinearForm(fes)
f += 32 * (y*(1-y)+x*(1-x)) * v * dx
a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx
a.Assemble()
f.Assemble()
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
Draw (gfu)


# POSTPROCESSING AND EXPORTING

# empty



# CLOSING STEPS

# empty



