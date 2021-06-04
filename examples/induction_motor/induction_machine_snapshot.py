from adze_modeler.geometry import Geometry
from adze_modeler.femm_wrapper import FemmWriter


def create_snapshot():
    """Create a snapshot manually, where every quantity can be calculated for one given position of the motor."""

    # the geometry is defined in a dxf file
    geo = Geometry()
    geo.import_dxf("2horse.dxf")
    print(geo)

    # create and save a model geometry in lua language
    femm_model = FemmWriter()
    femm_model.lua_model = femm_model.init_problem()
    femm_model.lua_model += femm_model.create_geometry(geo)
    femm_model.write('2horse.lua')
    return


if __name__ == "__main__":
    create_snapshot()
