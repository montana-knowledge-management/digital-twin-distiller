import fenics as fn
from adze_modeler.geometry import Geometry
from adze_modeler.gmsh import GMSHModel
from adze_modeler.boundaries import BoundaryCondition

from importlib_resources import files
import pyvista as pv
import meshio

from dolfin import Mesh, XDMFFile, File, MeshValueCollection, cpp, Measure,\
                   DirichletBC, FunctionSpace, Constant, TrialFunction, \
                   TestFunction, dot, grad, Function, solve, plot, dx

import matplotlib.pyplot as plt

def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read": [cell_data]})
    if prune_z:
        out_mesh.prune_z_0()
    return out_mesh


# importing the hand-made svg-image to use it as a FEM model
eml = files("examples.fenics_poisson_subdomain").joinpath("poisson_domain.svg")
geo = Geometry()
geo.import_svg(eml.as_posix())

# set the tolerance to merge the given lines
geo.epsilon = 1e-6
geo.merge_points()
gmsh = GMSHModel(geo)

GND = '1'
V0 = '2'
# gnd, v0
gnd_dirichlet = [(47.6, -200.0), (97.6, -150.0), (100, -250.0), (150, -165.0), (148, -220.0)]
v0_dirichlet = [(87.0, -180.0), (120.0, -175.0), (120.0, -194.0)]

# assign boundary for the <gnd> description
lines = {geo.lines[i].id: geo.lines[i] for i in range(0, len(geo.lines))}

# bounds - dolphin cannot handle strings only numbers
gmsh.boundaries[GND] = []
gmsh.boundaries[V0] = []

for label in gnd_dirichlet:
    closest_line = min(lines.values(), key=lambda li: li.distance_to_point(label[0], label[1]))
    gmsh.boundaries[GND].append(closest_line.id)

# assign the potential for <v0> description
for label in v0_dirichlet:
    closest_line = min(lines.values(), key=lambda li: li.distance_to_point(label[0], label[1]))
    gmsh.boundaries[V0].append(closest_line.id)

print(gmsh.boundaries)

gmsh.lcar = 15  # the characteristic length of the applied mesh can be set manually in the gmsh class
gmsh.gmsh_writer('poisson_domain')
print('test', gmsh.boundary_queue_gmsh)

## plotting out the resulting mesh by pyvista
msh = pv.read('poisson_domain.msh')
msh.plot(show_edges=True, cpos="xy")

epsilon_0 = 8.85e-12

# partially working version
# # create a mesh for dolphin
# mesh_from_file = meshio.read("poisson_domain.msh")
# line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
# meshio.write("poisson_facet.xdmf", line_mesh)
#
# cell_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
# meshio.write("poisson_domain.xml", cell_mesh)
#
# mesh = Mesh("poisson_domain.xml")
# xdmf = XDMFFile(mesh.mpi_comm(), "poisson_domain.xml")
# xdmf.write(mesh)
# plot(mesh)


#### Convert files into gmsh mesh file format
Wri_path = './Dolfin/'

mesh_from_file = meshio.read("poisson_domain.msh")

cell_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
meshio.write(Wri_path+"mesh.xdmf", cell_mesh)

line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
meshio.write(Wri_path+"mf.xdmf", line_mesh)

#### import the dolphin mesh and transform again
# write out the mesh
mesh = fn.Mesh()
mvc = fn.MeshValueCollection("size_t", mesh, 2) # 2 should change with fulldim parameter
with fn.XDMFFile(Wri_path+"mesh.xdmf") as infile:
    infile.read(mesh)
    infile.read(mvc, "name_to_read")
File(Wri_path+"Dolfin_circle_mesh.pvd").write(mesh)
mf = fn.cpp.mesh.MeshFunctionSizet(mesh, mvc)

### import boundaries -- fulldim - 1
mvc2 = fn.MeshValueCollection("size_t", mesh, 1)
with fn.XDMFFile(Wri_path+"mf.xdmf") as infile:
    infile.read(mvc2, "name_to_read")
cf = fn.cpp.mesh.MeshFunctionSizet(mesh, mvc2)

#File(Wri_path+"Dolfin_circle_facets.pvd").write(mvc2)

# plot(mesh)
# # Hold plot
# plt.show()

### FEM solution ###
V = FunctionSpace(mesh, 'Lagrange', 1)

# Define boundary conditions base on GMSH mesh marks [Physical Curves: 1, 2, 3, 4]
bc1 = fn.DirichletBC(V, fn.Constant(10.0), cf, 3)
bc2 = DirichletBC(V, Constant(0.0), cf, 2)

bc = [bc1, bc2]

# Define variational problem
u = fn.TrialFunction(V)
v = fn.TestFunction(V)

a = fn.dot(fn.grad(u), fn.grad(v)) * fn.dx
L = fn.Constant('0') * v * fn.dx
u = fn.Function(V)
fn.solve(a == L, u, bc)

File("Solution.pvd").write(u)

electric_field = fn.project(-fn.grad(u))
# solution

plt.subplot(1,2,1)
fn.plot(mf, 'domains')

plt.subplot(1,2,2)
fn.plot(u, 'Potential')
fn.plot(electric_field, title='Electric field')

plt.show()

#### final plots with pyvista

#fm = pv.read("Solution.pvd")
#qual = fm.compute_cell_quality()
#fm.plot(show_edges=True, cpos="xy")

