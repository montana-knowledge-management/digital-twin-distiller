from adze_modeler.objects import Node
from adze_modeler.femm_wrapper import FemmWriter, femm_fields


# the label postions  will be stored in the material class
# more then one region can be assigned to a material

class MaterialSnapshot:
    """
    The goal of this module is to insert the material label definitions according to the selected FEM software.
    This snapshot can be exported from the material editor class or it can be defined separately.
    It contains a region label list, which generates the labels to the defined regions.
    """

    def __init__(self):
        self.field_type = None
        self.meshsize = None  # if not defined, it tries to create an automesh
        self.material_definition = None
        self.region_labels = []  # contains nodes, for every defined regions

    def create_block_label(self):
        """Create the material definition labels for every given region. """
        cmds = []

        # in case of a FEMM model
        cmds.append(FemmWriter().add_material(self.material_definition))

        if self.field_type in femm_fields:
            for point in self.region_labels:
                cmds.append(FemmWriter().add_blocklabel(point.x, point.y))
                cmds.append(FemmWriter().select_label(point.x, point.y))

                #if self.meshsize:
                #    automesh = 0
                #else:
                #    automesh = 1

                #cmds.append(FemmWriter().set_blockprop(blockname=self.material_definition.material_name,
                #                                       automesh=automesh,
                #                                       meshsize=self.meshsize, group=0))
                cmds.append(FemmWriter().clear_selected())
        return cmds
