from adze_modeler.agros2d_wrapper import Agros2DWrapper

# mm
sc = 1e-3
Da = 31 * sc
Df = 2.5 * sc
Dc = 4 * sc
Dp = 15 * sc
Dg = 0.3 * sc
Hft = 2 * sc
Hg = 3 * sc
Hpm = 1 * sc
Hm = 22 * sc
Hf = 30 * sc
Wbox = 240 * sc
Hbox = 120 * sc

ag = Agros2DWrapper()
ag.add_field("m")
ag.set_coordinate_type("axisymmetric")
ag.field.set_solver("newton")

ag.field.add_boundary_condition("d", "a0", 0)

ag.field.add_material("air")
ag.field.add_material("aluminium")
# TODO: Implement nonlinear magnetic material
# ag.field.add_material("steel", mur="1;piecewise_linear,1,1$0.1,0.3,0.5;100,447,640")
ag.field.add_material("steel", mur=2364)
ag.field.add_material("Alnico5", remanence=1.28, remanence_angle=90)

ag.field.add_material("Coil_Upper", total_current=0)
ag.field.add_material("Coil_Lower", total_current=0)

ag.add_edge(Dg / 2, Hf / 2, Da / 2, Hf / 2)
ag.add_edge(Dg / 2, Hf / 2 - Hft, Da / 2 - Df - Dc, Hf / 2 - Hft)
ag.add_edge(Da / 2 - Df - Dc, Hf / 2 - Hft, Da / 2 - Df, Hf / 2 - Hft)
ag.add_edge(Da / 2 - Df, Hf / 2 - Hft, Da / 2, Hf / 2 - Hft)
ag.add_edge(0, Hm / 2, Dp / 2, Hm / 2)
ag.add_edge(Da / 2 - Df - Dc, Hpm / 2, Da / 2 - Df, Hpm / 2)

ag.add_edge(Dg / 2, -Hf / 2, Da / 2, -Hf / 2)
ag.add_edge(Dg / 2, -Hf / 2 + Hft, Da / 2 - Df - Dc, -Hf / 2 + Hft)
ag.add_edge(Da / 2 - Df - Dc, -Hf / 2 + Hft, Da / 2 - Df, -Hf / 2 + Hft)
ag.add_edge(Da / 2 - Df, -Hf / 2 + Hft, Da / 2, -Hf / 2 + Hft)
ag.add_edge(0, -Hm / 2, Dp / 2, -Hm / 2)
ag.add_edge(Da / 2 - Df - Dc, -Hpm / 2, Da / 2 - Df, -Hpm / 2)


ag.add_edge(Dg / 2, Hf / 2, Dg / 2, Hf / 2 - Hft)
ag.add_edge(Dg / 2, -Hf / 2, Dg / 2, -Hf / 2 + Hft)
ag.add_edge(0, Hm / 2, 0, -Hm / 2, boundary="a0")
ag.add_edge(Dp / 2, Hm / 2, Dp / 2, -Hm / 2)
ag.add_edge(Da / 2 - Df - Dc, Hf / 2 - Hft, Da / 2 - Df - Dc, Hpm / 2)
ag.add_edge(Da / 2 - Df, Hf / 2 - Hft, Da / 2 - Df, Hpm / 2)
ag.add_edge(Da / 2 - Df - Dc, -Hf / 2 + Hft, Da / 2 - Df - Dc, -Hpm / 2)
ag.add_edge(Da / 2 - Df, -Hf / 2 + Hft, Da / 2 - Df, -Hpm / 2)
ag.add_edge(Da / 2, Hf / 2, Da / 2, Hf / 2 - Hft)
ag.add_edge(Da / 2, -Hf / 2, Da / 2, -Hf / 2 + Hft)
ag.add_edge(Da / 2, Hf / 2 - Hft, Da / 2, -Hf / 2 + Hft)
ag.add_edge(Da / 2 - Df - Dc, Hpm / 2, Da / 2 - Df - Dc, -Hpm / 2)
ag.add_edge(Da / 2 - Df, Hpm / 2, Da / 2 - Df, -Hpm / 2)

ag.add_edge(0, Hm / 2, 0, Hbox / 2, boundary="a0")
ag.add_edge(0, -Hm / 2, 0, -Hbox / 2, boundary="a0")
ag.add_edge(0, -Hbox / 2, Wbox / 2, -Hbox / 2, boundary="a0")
ag.add_edge(0, Hbox / 2, Wbox / 2, Hbox / 2, boundary="a0")
ag.add_edge(Wbox / 2, Hbox / 2, Wbox / 2, -Hbox / 2, boundary="a0")


ag.add_block_label(Dp / 4, 0, "Alnico5")
ag.add_block_label(Da / 2 - Df / 2, 0, "steel")
ag.add_block_label(Wbox / 4, 0, "air")
ag.add_block_label(Da / 2 - Df - Dc / 2, 0, "aluminium")
ag.add_block_label(Da / 4, Hf / 2 - Hft / 2, "aluminium")
ag.add_block_label(Da / 4, -Hf / 2 + Hft / 2, "aluminium")
ag.add_block_label(Da / 2 - Df - Dc / 2, (Hf / 2 - Hft) / 2, "Coil_Upper")
ag.add_block_label(Da / 2 - Df - Dc / 2, -(Hf / 2 - Hft) / 2, "Coil_Lower")

ag.export("a2g_magnetic.py")

from subprocess import run

run(["agros2d", "-s", "a2g_magnetic.py"])
