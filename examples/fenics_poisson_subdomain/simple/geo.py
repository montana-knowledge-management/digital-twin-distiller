import pygmsh
resolution = 0.01

# Channel parameters
L = 2.2
H = 0.41
c = [0.2, 0.2, 0]
r = 0.05

# Initialize empty geometry using the build in kernel in GMSH
geometry = pygmsh.geo.Geometry()
# Fetch model we would like to add data to
model = geometry.__enter__()
# Add circle
circle = model.add_circle(c, r, mesh_size=resolution)

# Add points with finer resolution on left side
points = [model.add_point((0, 0, 0), mesh_size=resolution),
          model.add_point((L, 0, 0), mesh_size=5*resolution),
          model.add_point((L, H, 0), mesh_size=5*resolution),
          model.add_point((0, H, 0), mesh_size=resolution)]

# Add lines between all points creating the rectangle
channel_lines = [model.add_line(points[i], points[i+1])
                 for i in range(-1, len(points)-1)]

# Create a line loop and plane surface for meshing
channel_loop = model.add_curve_loop(channel_lines)
plane_surface = model.add_plane_surface(
    channel_loop, holes=[circle.curve_loop])

# Call gmsh kernel before add physical entities
model.synchronize()

volume_marker = 6
model.add_physical([plane_surface], "Volume")
model.add_physical([channel_lines[0]], "Inflow")
model.add_physical([channel_lines[2]], "Outflow")
model.add_physical([channel_lines[1], channel_lines[3]], "Walls")
model.add_physical(circle.curve_loop.curves, "Obstacle")


# We generate the mesh using the pygmsh function `generate_mesh`. Generate mesh returns a `meshio.Mesh`. However, this mesh is tricky to extract physical tags from. Therefore we write the mesh to file using the `gmsh.write` function.
geometry.generate_mesh(dim=2)
import gmsh
gmsh.write("mesh.msh")
gmsh.clear()
geometry.__exit__()


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