from adze_modeler.geometry import Geometry
from adze_modeler.gmsh import GMSHModel
from adze_modeler.objects import Node, Line

import pyvista as pv


# define the geometry by hand, a simple arc
geo = Geometry()

a = Node(x=0.0, y=0.0, id=1)
b = Node(x=10.0, y=0.0, id=2)
c = Node(x=0.0, y=10.0, id=3)

geo.add_line(Line(a, b))
geo.add_line(Line(a, c))
geo.add_line(Line(c, b))

# adze can find the surface automatically by the built-in commands
surfaces = geo.find_surfaces()

# the connection graph of the connecting nodes:
geo.plot_connection_graph()

# the gmsh-writer class can automatically save the given graph into vtk format
gmsh = GMSHModel(geo)
gmsh.lcar = 0.1  # the characteristic length of the applied mesh can be set manually in the gmsh class
gmsh.gmsh_writer('triangle')

# plotting out the resulting mesh by pyvista
msh = pv.read('triangle.vtk')
msh.plot(show_edges=True, cpos="xy")
