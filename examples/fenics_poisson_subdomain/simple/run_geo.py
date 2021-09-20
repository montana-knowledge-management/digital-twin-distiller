# start
print('start')

# import
import fenics as fn
import dolfin
import sys
import numpy
import matplotlib.pyplot as plt
import pygmsh
import gmsh

fullDim = 2 #2D model

mesh = fn.Mesh()
mvc = fn.MeshValueCollection("size_t", mesh, fullDim)
with fn.XDMFFile("mesh.xdmf") as infile:
   infile.read(mesh)
   infile.read(mvc, "name_to_read")
mf = fn.cpp.mesh.MeshFunctionSizet(mesh, mvc)

mvc2 = fn.MeshValueCollection("size_t", mesh, fullDim-1)
with fn.XDMFFile("facet_mesh.xdmf") as infile:
   infile.read(mvc2, "name_to_read")
cf = fn.cpp.mesh.MeshFunctionSizet(mesh, mvc2)


# Define function space
V = fn.FunctionSpace(mesh, 'Lagrange', 1)
#markers = fn.MeshFunction('size_t', mesh, fullDim, cf)
#dx = fn.Measure('dx', domain=mesh, subdomain_data=markers)
#dx = fn.Measure('dx', domain=mesh, subdomain_data=cf)
dx = fn.Measure('dx', domain=mesh, subdomain_data=mf)

# Define boundary condition
bc = fn.DirichletBC(V, fn.Constant(0), 'on_boundary')
bcs = [bc]
bc = fn.DirichletBC(V, fn.Constant(1), cf, 5)
bcs.append(bc)

# some plot setup
plt.figure()
plt.subplot(1,3,1)
fn.plot(mf, title='Sub Domains')


# equation solving
u = fn.TrialFunction(V)
v = fn.TestFunction(V)
a = fn.dot(fn.grad(u), fn.grad(v)) * fn.dx
L = fn.Constant('0') * v * fn.dx
u = fn.Function(V)
fn.solve(a == L, u, bcs)

electric_field = fn.project(-fn.grad(u))
plt.subplot(1,3,2)
fn.plot(u, title='Fields')
plt.subplot(1,3,3)
fn.plot(electric_field, title='Solution')

plt.show()