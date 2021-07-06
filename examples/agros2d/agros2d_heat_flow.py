from adze_modeler.agros2d_wrapper import Agros2DWrapper

ag = Agros2DWrapper()
ag.add_field("h")


ag.export("a2g_heat_flow.py")
