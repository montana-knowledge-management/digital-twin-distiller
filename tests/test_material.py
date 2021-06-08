import unittest

from adze_modeler.femm_wrapper import femm_magnetic
from adze_modeler.femm_wrapper import MagneticMaterial
from adze_modeler.material import MaterialSnapshot
from adze_modeler.objects import Node


class TestMaterialSnapshot(unittest.TestCase):
    def test_add_femm_magnetic_material(self):
        coil = MagneticMaterial("coil", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)

        # creates a snapshot for 3 coil regions in the geometry with 'coil' material, automatic meshsize
        coils = MaterialSnapshot()
        coils.material_definition = coil
        coils.field_type = femm_magnetic  # magnetic field with the femm material

        # selected regions
        a = Node(x=0.5, y=1.00, id=15)
        b = Node(x=5.0, y=10.0, id=16)
        c = Node(x=2.5, y=3.0, id=17)

        coils.region_labels = [a, b, c]

        # generates the essential lua script commands
        cmd = coils.create_femm_block_label()

        self.assertEqual(cmd[1], "mi_addblocklabel(0.5, 1.0)")
        self.assertEqual(cmd[2], "mi_selectlabel(0.5, 1.0)")
        self.assertEqual(cmd[0], "mi_addmaterial('coil', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)")
        # self.assertEqual(cmd[3], "mi_setblockprop('coil', 1, None, 'incircuit', 0, 0, 0)")
