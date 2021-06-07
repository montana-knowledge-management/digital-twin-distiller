from adze_modeler.geometry import Geometry
from adze_modeler.objects import Node
from adze_modeler.femm_wrapper import FemmWriter, MagneticMaterial, femm_magnetic
from adze_modeler.material import MaterialSnapshot

from math import pi


def create_snapshot():
    """Create a snapshot manually, where every quantity can be calculated for one given position of the motor."""

    # the geometry is defined in a dxf file
    geo = Geometry()
    geo.import_dxf("2horse.dxf")

    # create and save a model geometry in lua language
    femm_model = FemmWriter()
    femm_model.lua_model = femm_model.init_problem()
    femm_model.lua_model += femm_model.create_geometry(geo)

    # material properties
    # silicon core with linea bh relationships
    silicon_core = MagneticMaterial("silicon_core", 7000, 7000, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
    aluminum = MagneticMaterial("aluminum", 1, 1, 0, 0, 34.5, 0, 0, 1, 0, 0, 0, 0, 0)
    coil = MagneticMaterial("copper", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
    air = MagneticMaterial("air", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)

    # airgap ----------------------------------------------------------------------
    airgap = MaterialSnapshot()
    airgap.material_definition = air
    airgap.field_type = femm_magnetic  # magnetic field with the femm material

    # label positions
    airgap.region_labels = [Node(40.5, 3.5)]
    femm_model.lua_model += airgap.create_block_label()

    # cores -----------------------------------------------------------------------
    cores = MaterialSnapshot()
    cores.material_definition = silicon_core
    cores.field_type = femm_magnetic  # magnetic field with the femm material

    # label positions
    cores.region_labels = [Node(16.0, 16.0), Node(40.0, 40.0)]
    femm_model.lua_model += cores.create_block_label()

    # aluminum rods ----------------------------------------------------------------
    rods = MaterialSnapshot()
    rods.material_definition = aluminum
    rods.field_type = femm_magnetic  # magnetic field with the femm material

    base_node = Node(31.6, 0.33)
    for i in range(0, 8):
        rods.region_labels.append(base_node)
        base_node = base_node.rotate(rad=12.5 / 180 * pi)
    femm_model.lua_model += rods.create_block_label()

    # save to file
    femm_model.write('2horse.lua')
    return


if __name__ == "__main__":
    create_snapshot()
