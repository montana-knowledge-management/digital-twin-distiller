import unittest

from adze_modeler.femm_wrapper import FemmWriter
from adze_modeler.femm_wrapper import FemmExecutor

# class TestFemmHeatProblem(unittest.TestCase):
#
#
#     def test_heat_problem(self):
#         writer = FemmWriter()
#         # TODO: FemmWriter.set_field(fieldtype) ?
#         writer.field = "heat_flow"
#         writer.lua_model.extend(writer.init_problem())
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
#         writer.lua_model.append(writer.add_node(-w2 / 2, -h2 / 2))
#         writer.lua_model.append(writer.add_node(w2 / 2, -h2 / 2))
#         writer.lua_model.append(writer.add_node(w2 / 2, h2 / 2))
#         writer.lua_model.append(writer.add_node(-w2 / 2, h2 / 2))
#
#         writer.lua_model.append(writer.add_segment(-w1 / 2, -h1 / 2, -w2 / 2, -h2 / 2))
#         writer.lua_model.append(writer.add_segment(-w2 / 2, -h2 / 2, -w2 / 2, h2 / 2))
#         writer.lua_model.append(writer.add_segment(-w2 / 2, h2 / 2, w2 / 2, h2 / 2))
#         writer.lua_model.append(writer.add_segment(w2 / 2, h2 / 2, w2 / 2, -h2 / 2))
#
#         writer.lua_model.extend(writer.close())
#         writer.write("test.lua")
#         FemmExecutor().run_femm("test.lua")
