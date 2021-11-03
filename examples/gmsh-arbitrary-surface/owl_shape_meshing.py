import pyvista as pv
from importlib_resources import files

from adze_modeler.geometry import Geometry
from adze_modeler.gmsh import GMSHModel

# importing the hand-made svg-image to use it as a FEM model
eml = files("examples.gmsh-arbitrary-surface.resources").joinpath(
    "owl-shape.svg"
)
geo = Geometry()
geo.import_svg(eml.as_posix())

# set the tolerance to merge the given lines
geo.epsilon = 1e-6
geo.merge_points()

# find_surfaces function can automatically detect the closed loops in the created geometry
# and find the user defined surfaces
surfaces = geo.find_surfaces()

# the connection graph of the connecting nodes:
geo.plot_connection_graph()

# the gmsh-writer class can automatically save the given graph into vtk format
gmsh = GMSHModel(geo)
gmsh.lcar = 10.0  # the characteristic length of the applied mesh can be set manually in the gmsh class
gmsh.gmsh_writer("owl_shape")

# plotting out the resulting mesh by pyvista
msh = pv.read("owl_shape.msh")
msh.plot(show_edges=True, cpos="xy")
