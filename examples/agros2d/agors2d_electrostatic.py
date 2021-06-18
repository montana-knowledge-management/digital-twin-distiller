from adze_modeler.agros2d_wrapper import Agros2DWrapper

a = 10  # cm
epsilon_r = 2.1  # -
Ug = 10  # V
# the length of the excitation
a1 = 0.1 # cm

ag = Agros2DWrapper()
ag.add_field('e')

ag.field.add_boundary_condition('d', 'Vg', Ug)
ag.field.add_boundary_condition('d', 'GND', 0)
ag.field.add_boundary_condition('n', 'neumann', 0)

ag.field.add_material("PVC", 4.0)
ag.add_block_label(a/4, a/4, "PVC")

ag.add_edge(-a/2, 0, -a1/2, 0, boundary="neumann")
ag.add_edge(-a1/2, 0, 0, 0, boundary="Vg")
ag.add_edge(0, 0, 0, -a1/2, boundary="Vg")
ag.add_edge(0, -a1/2, 0, -a/2, boundary="neumann")
ag.add_edge(0, -a/2, a/2, -a/2, boundary="neumann")
ag.add_edge(a / 2, -a/2, a/2, a/2, boundary="GND")
ag.add_edge(a / 2, a/2, -a/2, a/2, boundary="GND")
ag.add_edge(-a / 2, a/2, -a/2, 0, boundary="neumann")

ag.export("a2g_electrostatic.py")