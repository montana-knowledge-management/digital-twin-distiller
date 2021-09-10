from adze_modeler.geometry import Geometry
from adze_modeler.gmsh import GMSHModel
from adze_modeler.boundaries import BoundaryCondition

from importlib_resources import files
import pyvista as pv

# importing the hand-made svg-image to use it as a FEM model
eml = files("examples.fenics_poisson_subdomain").joinpath("poisson_domain.svg")
geo = Geometry()
geo.import_svg(eml.as_posix())

# set the tolerance to merge the given lines
geo.epsilon = 1e-6
geo.merge_points()
gmsh = GMSHModel(geo)

gmsh.label_queue.append((47.6, -200.0, "gnd"))
gmsh.label_queue.append((87.0, -180.0, "v0"))
# assign boundary for the <gnd> description
lines = {geo.lines[i].id: geo.lines[i] for i in range(0, len(geo.lines))}
closest_line = min(lines.values(), key=lambda li: li.distance_to_point(47.6, -200.0))
gmsh.boundaries['gnd'] = [closest_line.id]

# assign the potential for <v0> description
closest_line = min(lines.values(), key=lambda li: li.distance_to_point(87.0, -180.0))
gmsh.boundaries['v0'] = [closest_line.id]
print(gmsh.boundaries)

gmsh.lcar = 10.0  # the characteristic length of the applied mesh can be set manually in the gmsh class
gmsh.gmsh_writer('poisson_domain')
print('test', gmsh.boundary_queue_gmsh)

# plotting out the resulting mesh by pyvista
msh = pv.read('poisson_domain.msh')
msh.plot(show_edges=True, cpos="xy")
