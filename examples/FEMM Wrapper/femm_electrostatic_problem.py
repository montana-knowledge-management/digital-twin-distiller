# import os
# import unittest
#
# from adze_modeler.femm_wrapper import ElectrostaticMaterial
# from adze_modeler.femm_wrapper import FemmExecutor
# from adze_modeler.femm_wrapper import FemmWriter
# from adze_modeler.femm_wrapper import kw_electrostatic
#
#
# class TestFemmElectrostaticProblem(unittest.TestCase):
#     # integration test ignored from the unittest list
#     def test_electrostatic_problem(self):
#         writer = FemmWriter()
#         writer.field = kw_electrostatic
#         writer.lua_model.extend(writer.init_problem("electrostatic_data.csv"))
#
#         writer.lua_model.append(writer.electrostatic_problem("centimeters", "planar"))
#
#         a = 10  # cm
#         epsilon_r = 2.1  # -
#         Ug = 10  # V
#
#         writer.lua_model.append(writer.add_node(-a / 2, a / 2))
#         writer.lua_model.append(writer.add_node(-a / 2, 0))
#         writer.lua_model.append(writer.add_node(0, 0))
#         writer.lua_model.append(writer.add_node(0, -a / 2))
#         writer.lua_model.append(writer.add_node(a / 2, -a / 2))
#         writer.lua_model.append(writer.add_node(a / 2, a / 2))
#
#         writer.lua_model.append(writer.add_segment(-a / 2, a / 2, -a / 2, 0))
#         writer.lua_model.append(writer.add_segment(-a / 2, 0, 0, 0))
#         writer.lua_model.append(writer.add_segment(0, 0, 0, -a / 2))
#         writer.lua_model.append(writer.add_segment(0, -a / 2, a / 2, -a / 2))
#         writer.lua_model.append(writer.add_segment(a / 2, -a / 2, a / 2, a / 2))
#         writer.lua_model.append(writer.add_segment(a / 2, a / 2, -a / 2, a / 2))
#
#         # Adding material properties
#         blocklabel = (a / 4, a / 4)
#         mat = ElectrostaticMaterial("Teflon", epsilon_r, epsilon_r, 0)
#         writer.lua_model.append(writer.add_material(mat))
#         writer.lua_model.append(writer.add_blocklabel(*blocklabel))
#         writer.lua_model.append(writer.select_label(*blocklabel))
#         writer.lua_model.append(writer.set_blockprop("Teflon"))
#
#         # Adding boundary properties
#         writer.lua_model.append(writer.add_pointprop("Ug", Vp=Ug))
#         writer.lua_model.append(writer.add_pointprop("U0", Vp=0))
#
#         writer.lua_model.append(writer.select_node(0, 0))
#         writer.lua_model.append(writer.set_pointprop("Ug"))
#         writer.lua_model.append("ei_clearselected()")
#
#         writer.lua_model.append(writer.select_node(a / 2, a / 2))
#         writer.lua_model.append(writer.set_pointprop("U0"))
#         writer.lua_model.append("ei_clearselected()")
#
#         writer.lua_model.append("ei_zoomnatural()")
#         writer.lua_model.append("ei_zoomout()")
#         writer.lua_model.append("hideconsole()")
#         writer.lua_model.append(writer.save_as("electrostatic_test.fee"))
#         writer.lua_model.append(writer.analyze())
#         writer.lua_model.append(writer.load_solution())
#
#         # Examine the results
#         writer.lua_model.append(f"eo_selectblock({a / 4}, {a / 4})")
#         writer.lua_model.append("E = eo_blockintegral(0)")  # Stored Energy
#         writer.lua_model.append(writer.write_out_result("E", "E"))
#         writer.lua_model.extend(writer.close())
#
#         # Execute the script
#         writer.write("electrostatic_test.lua")
#         FemmExecutor().run_femm("../integration_tests/electrostatic_test.lua")
#
#         with open("electrostatic_data.csv") as f:
#             content = f.readlines()
#             E = content[0].split(",")[1]  # Stored energy
#             self.assertEqual(round(float(E), 4), round(1.5293 * 10 ** -12, 4))
#
#         # try:
#         #     os.remove("electrostatic_data.csv")
#         #     os.remove("electrostatic_test.fee")
#         #     os.remove("electrostatic_test.lua")
#         #     os.remove("electrostatic_test.res")
#         # except FileNotFoundError:
#         #     pass
#
#
# if __name__ == '__main__':
#     TestFemmElectrostaticProblem.integation_test_electrostatic_problem()
