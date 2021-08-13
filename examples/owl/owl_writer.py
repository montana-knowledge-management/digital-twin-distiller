from adze_modeler.geometry import Geometry
from adze_modeler.gmsh import GMSHModel

from importlib_resources import files

# open the owl svg
eml = files("examples.owl").joinpath("owl-shape.svg")
geo = Geometry()
geo.import_svg(eml.as_posix())
# set the tolerance to merge the given lines
geo.epsilon = 1e-6
print(geo)
#geo.merge_points()

geo.find_surfaces()
geo.export_svg('bagoly.svg')
# create a gmsh mesh from the given geometry
gmsh = GMSHModel(geo)
#gmsh.add_lines()
#gmsh.gmsh_geometry.save_geometry("1.geo_unrolled")
gmsh.gmsh_writer()