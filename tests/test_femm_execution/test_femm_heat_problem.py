# import os
# import unittest
#
# from adze_modeler.femm_wrapper import FemmExecutor
# from adze_modeler.femm_wrapper import FemmWriter
# from adze_modeler.femm_wrapper import HeatFlowConvection
# from adze_modeler.femm_wrapper import HeatFlowFixedTemperature
# from adze_modeler.femm_wrapper import HeatFlowMaterial
# from adze_modeler.femm_wrapper import kw_heat_flow
#
#
# def c2k(C):
#     return C + 274.15
#
#
# class TestFemmHeatProblem(unittest.TestCase):
#     def test_heat_problem(self):
#         writer = FemmWriter()
#         writer.field = kw_heat_flow
#         writer.lua_model.extend(writer.init_problem("heat_data.csv"))
#
#         writer.lua_model.append(writer.heat_problem("meters", "planar"))
#
#         w1 = 0.3 + 0.2 + 0.3
#         h1 = 0.4
#         w2 = 0.2
#         h2 = 0.2
#
#         # Constructing the geometry
#         writer.lua_model.append(writer.add_node(-w1 / 2, -h1 / 2))
#         writer.lua_model.append(writer.add_node(w1 / 2, -h1 / 2))
#         writer.lua_model.append(writer.add_node(w1 / 2, h1 / 2))
#         writer.lua_model.append(writer.add_node(-w1 / 2, h1 / 2))
#
#         writer.lua_model.append(writer.add_node(-w2 / 2, -h1 / 2))
#         writer.lua_model.append(writer.add_node(-w2 / 2, -h1 / 2 + h2))
#         writer.lua_model.append(writer.add_node(w2 / 2, -h1 / 2 + h2))
#         writer.lua_model.append(writer.add_node(w2 / 2, -h1 / 2))
#
#         # adding segments
#         writer.lua_model.append(writer.add_segment(-w1 / 2, -h1 / 2, -w2 / 2, -h1 / 2))
#         writer.lua_model.append(writer.add_segment(-w2 / 2, -h1 / 2, -w2 / 2, -h1 / 2 + h2))
#         writer.lua_model.append(writer.add_segment(-w2 / 2, -h1 / 2 + h2, w2 / 2, -h1 / 2 + h2))
#         writer.lua_model.append(writer.add_segment(w2 / 2, -h1 / 2 + h2, w2 / 2, -h1 / 2))
#         writer.lua_model.append(writer.add_segment(w2 / 2, -h1 / 2, w1 / 2, -h1 / 2))
#         writer.lua_model.append(writer.add_segment(w1 / 2, -h1 / 2, w1 / 2, h1 / 2))
#         writer.lua_model.append(writer.add_segment(w1 / 2, h1 / 2, -w1 / 2, h1 / 2))
#         writer.lua_model.append(writer.add_segment(-w1 / 2, h1 / 2, -w1 / 2, -h1 / 2))
#
#         # Adding material properties
#         blocklabel = (0.0, h1 / 2.0 * 0.6)
#         mat = HeatFlowMaterial("Material", 1380, 1380, 0, 0)
#         writer.lua_model.append(writer.add_material(mat))
#         writer.lua_model.append(writer.add_blocklabel(*blocklabel))
#         writer.lua_model.append(writer.select_label(*blocklabel))
#         writer.lua_model.append(writer.set_blockprop("Material"))
#
#         # Adding boundary properties
#         insulation = HeatFlowConvection("ins", 0, c2k(22))
#         convection = HeatFlowConvection("conv", 50, c2k(50))
#         fixtemp = HeatFlowFixedTemperature("fixtemp", c2k(300))
#         writer.lua_model.append(writer.add_boundary(insulation))
#         writer.lua_model.append(writer.add_boundary(convection))
#         writer.lua_model.append(writer.add_boundary(fixtemp))
#
#         # Set boundary properties
#
#         # insulation
#         writer.lua_model.append(writer.select_segment(-w1 / 2, 0))
#         writer.lua_model.append(writer.select_segment(w1 / 2, 0))
#         writer.lua_model.append(writer.set_segment_prop("ins"))
#         writer.lua_model.append("hi_clearselected()")
#
#         # convection
#         writer.lua_model.append(writer.select_segment(0, h1 / 2))
#         writer.lua_model.append(writer.set_segment_prop("conv"))
#         writer.lua_model.append("hi_clearselected()")
#
#         # fixed temperature
#         writer.lua_model.append(writer.select_segment(-(w2 - w1) / 4, -h1 / 2))
#         writer.lua_model.append(writer.select_segment((w2 - w1) / 4, -h1 / 2))
#         writer.lua_model.append(writer.select_segment(-w2 / 2, -h1 / 2 + h2 / 2))
#         writer.lua_model.append(writer.select_segment(w2 / 2, -h1 / 2 + h2 / 2))
#         writer.lua_model.append(writer.select_segment(0, -h1 / 2 + h2))
#         writer.lua_model.append(writer.set_segment_prop("fixtemp"))
#         writer.lua_model.append("hi_clearselected()")
#
#         writer.lua_model.append("hi_zoomnatural()")
#         writer.lua_model.append("hi_zoomout()")
#         writer.lua_model.append("hideconsole()")
#         writer.lua_model.append(writer.save_as("heatflow_test.feh"))
#         writer.lua_model.append(writer.analyze())
#         writer.lua_model.append(writer.load_solution())
#
#         # select the material
#         writer.lua_model.append("ho_selectblock(0, 0.1)")
#         writer.lua_model.append("Fx, Fy = ho_blockintegral(3)")
#         writer.lua_model.append(writer.write_out_result("Fx", "Fx"))
#         writer.lua_model.append(writer.write_out_result("Fy", "Fy"))
#         writer.lua_model.extend(writer.close())
#
#         writer.write("heatflow_test.lua")
#         FemmExecutor().run_femm("heatflow_test.lua")
#
#         with open("heat_data.csv") as f:
#             content = f.readlines()
#             # print(content[0])
#             # print(content[1])
#             Fx = content[0].split(",")[1]
#             Fy = content[1].split(",")[1]
#             self.assertEqual(round(float(Fx), 4), 0.0112)  # 0.0112
#             self.assertEqual(round(float(Fy), 4), 9551.0549)  # 9551.05
#         try:
#             with open("heat_data.csv") as f:
#                 content = f.readlines()
#                 # print(content[0])
#                 # print(content[1])
#                 Fx = content[0].split(",")[1]
#                 Fy = content[1].split(",")[1]
#                 self.assertEqual(round(float(Fx), 4), 0.0112)  # 0.0112
#                 self.assertEqual(round(float(Fy), 4), 9551.0549)  # 9551.05
#
#             os.remove("heat_data.csv")
#             os.remove("heatflow_test.anh")
#             os.remove("heatflow_test.feh")
#             os.remove("heatflow_test.lua")
#
#         except FileNotFoundError:
#             self.assertTrue(False)
