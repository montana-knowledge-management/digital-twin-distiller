from adze_modeler.geometry import Geometry
from adze_modeler.gmsh import GMSHModel

from importlib_resources import files

# test - 0
eml = files("tests.pygmsh_tests").joinpath("test1.svg")
geo = Geometry()
geo.import_svg(eml.as_posix())
# set the tolerance to merge the given lines
geo.epsilon = 1e-6
geo.merge_points()
#print(geo)

surfaces = geo.find_surfaces()
#print('surfaces', surfaces)
geo.plot_connection_graph()
geo.export_svg('test1.svg')
# create a gmsh mesh from the given geometry
gmsh = GMSHModel(geo)
gmsh.gmsh_writer('test1')