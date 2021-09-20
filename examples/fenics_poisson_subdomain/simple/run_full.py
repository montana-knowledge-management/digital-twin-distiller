

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

##### this is new ################
import meshio
mesh_from_file = meshio.read("mesh.msh")

import numpy
def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
    if prune_z:
        out_mesh.prune_z_0()
    return out_mesh

line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
meshio.write("facet_mesh.xdmf", line_mesh)

triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
meshio.write("mesh.xdmf", triangle_mesh)

######### the rest is unchanged ###############

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

## and then the rest
