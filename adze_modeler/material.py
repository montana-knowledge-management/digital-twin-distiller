from adze_modeler.objects import Node
from adze_modeler.femm_wrapper import FemmWriter, femm_magnetic, femm_fields


# the label positions  will be stored in the material class
# more then one region can be assigned to a material

class MaterialSnapshot:
    """
    The goal of this module is to insert the material label definitions according to the selected FEM software.
    This snapshot can be exported from the material editor class or it can be defined separately.
    It contains a region label list, which generates the labels to the defined regions.
    """

    def __init__(self):

        # required fields
        self.field_type = None
        self.meshsize = None  # if not defined, it tries to create an automesh
        self.material_definition = None
        self.region_labels = []  # contains nodes, for every defined regions

        # optional parameters
        self.group = 0

        # optional femm parameters, which has to be defined only in the case of magnetic simulation
        self.circuit_name = "<None>"
        self.turn_umber = 0
        self.magdirection = 0
        self.circuit_type = 0  # parallel or series
        self.circuit_current = 0  # Amper

    def create_femm_block_label(self):
        """
            Creates the material definition labels for every given region.
            The circuit properties and material property definitions are merged within this commands.
        """
        cmds = []

        # in case of a FEMM model
        cmds.append(FemmWriter().add_material(self.material_definition))

        if self.field_type in femm_fields:
            for point in self.region_labels:
                cmds.append(FemmWriter().add_blocklabel(point.x, point.y))
                cmds.append(FemmWriter().select_label(point.x, point.y))

                if self.meshsize:
                    automesh = 0
                else:
                    automesh = 1

                # femm command contains special inputs only in the case of magnetic field simulation
                if self.field_type == femm_magnetic:
                    if self.circuit_name != "<None>":
                        cmds.append(FemmWriter().add_circprop(circuitname=self.circuit_name,
                                                              i=self.circuit_current,
                                                              circuittype=self.circuit_type))

                    cmds.append(FemmWriter().set_blockprop(blockname=self.material_definition.material_name,
                                                           automesh=automesh,
                                                           meshsize=self.meshsize,
                                                           group=0,
                                                           circuit_name=self.circuit_name,
                                                           magdirection=self.magdirection,
                                                           turns=self.turn_umber))
                else:
                    cmds.append(FemmWriter().set_blockprop(blockname=self.material_definition.material_name,
                                                           automesh=automesh,
                                                           meshsize=self.meshsize, group=0))
                cmds.append(FemmWriter().clear_selected())
        return cmds
