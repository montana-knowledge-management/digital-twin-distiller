import meshio
from dolfin import Mesh, XDMFFile, File, MeshValueCollection, cpp, Measure,\
                   DirichletBC, FunctionSpace, Constant, TrialFunction, \
                   TestFunction, dot, grad, dx, Function, solve

msh = meshio.read("./Square.msh")
Wri_path = './Dolfin_mesh_functions/'
meshio.write(Wri_path+"mesh.xdmf",
             meshio.Mesh(points=msh.points,
                         cells={"triangle": msh.cells["triangle"]}))

meshio.write(Wri_path+"mf.xdmf",
             meshio.Mesh(points=msh.points,
                         cells={"line": msh.cells["line"]},
                         cell_data={"line": {"name_to_read": msh.cell_data["line"]["gmsh:physical"]}}))

mesh = Mesh()
with XDMFFile(Wri_path+"mesh.xdmf") as infile:
    infile.read(mesh)
File(Wri_path+"Dolfin_circle_mesh.pvd").write(mesh)

mvc = MeshValueCollection("size_t", mesh, 1)
with XDMFFile(Wri_path+"mf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
File(Wri_path+"Dolfin_circle_facets.pvd").write(mf)

V = FunctionSpace(mesh, 'P', 1)

# Define boundary conditions base on GMSH mesh marks [Physical Curves: 1, 2, 3, 4]
bc1 = DirichletBC(V, Constant(0.0), mf, 1)
bc2 = DirichletBC(V, Constant(10.0), mf, 2)
bc3 = DirichletBC(V, Constant(6.0), mf, 3)
bc4 = DirichletBC(V, Constant(3.0), mf, 4)

bc = [bc1, bc2, bc3, bc4]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-15.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)
File("Solution.pvd").write(u)