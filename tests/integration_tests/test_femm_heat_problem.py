import unittest
from collections import Counter

from importlib_resources import files

from adze_modeler.femm_wrapper import (
    FemmWriter,
    HeatFlowConvection,
    HeatFlowFixedTemperature,
    HeatFlowMaterial,
    femm_heat_flow,
)


def c2k(C):
    return C + 274.15


class TestFemmCurrentFlowProblem(unittest.TestCase):
    def test_heat_problem(self):
        writer = FemmWriter()
        writer.field = femm_heat_flow
        writer.init_problem("heat_data.csv")

        writer.heat_problem("meters", "planar")

        w1 = 0.3 + 0.2 + 0.3
        h1 = 0.4
        w2 = 0.2
        h2 = 0.2

        # Constructing the geometry
        writer.add_node(-w1 / 2, -h1 / 2)
        writer.add_node(w1 / 2, -h1 / 2)
        writer.add_node(w1 / 2, h1 / 2)
        writer.add_node(-w1 / 2, h1 / 2)

        writer.add_node(-w2 / 2, -h1 / 2)
        writer.add_node(-w2 / 2, -h1 / 2 + h2)
        writer.add_node(w2 / 2, -h1 / 2 + h2)
        writer.add_node(w2 / 2, -h1 / 2)

        # adding segments
        writer.add_segment(-w1 / 2, -h1 / 2, -w2 / 2, -h1 / 2)
        writer.add_segment(-w2 / 2, -h1 / 2, -w2 / 2, -h1 / 2 + h2)
        writer.add_segment(-w2 / 2, -h1 / 2 + h2, w2 / 2, -h1 / 2 + h2)
        writer.add_segment(w2 / 2, -h1 / 2 + h2, w2 / 2, -h1 / 2)
        writer.add_segment(w2 / 2, -h1 / 2, w1 / 2, -h1 / 2)
        writer.add_segment(w1 / 2, -h1 / 2, w1 / 2, h1 / 2)
        writer.add_segment(w1 / 2, h1 / 2, -w1 / 2, h1 / 2)
        writer.add_segment(-w1 / 2, h1 / 2, -w1 / 2, -h1 / 2)

        # Adding material properties
        blocklabel = (0.0, h1 / 2.0 * 0.6)
        mat = HeatFlowMaterial("Material", 1380, 1380, 0, 0)
        writer.add_material(mat)
        writer.add_blocklabel(*blocklabel)
        writer.select_label(*blocklabel)
        writer.set_blockprop("Material")

        # Adding boundary properties
        insulation = HeatFlowConvection("ins", 0, c2k(22))
        convection = HeatFlowConvection("conv", 50, c2k(50))
        fixtemp = HeatFlowFixedTemperature("fixtemp", c2k(300))
        writer.add_boundary(insulation)
        writer.add_boundary(convection)
        writer.add_boundary(fixtemp)

        # Set boundary properties

        # insulation
        writer.select_segment(-w1 / 2, 0)
        writer.select_segment(w1 / 2, 0)
        writer.set_segment_prop("ins")
        writer.lua_model.append("hi_clearselected()")

        # convection
        writer.select_segment(0, h1 / 2)
        writer.set_segment_prop("conv")
        writer.lua_model.append("hi_clearselected()")

        # fixed temperature
        writer.select_segment(-(w2 - w1) / 4, -h1 / 2)
        writer.select_segment((w2 - w1) / 4, -h1 / 2)
        writer.select_segment(-w2 / 2, -h1 / 2 + h2 / 2)
        writer.select_segment(w2 / 2, -h1 / 2 + h2 / 2)
        writer.select_segment(0, -h1 / 2 + h2)
        writer.set_segment_prop("fixtemp")
        writer.lua_model.append("hi_clearselected()")

        writer.lua_model.append("hi_zoomnatural()")
        writer.lua_model.append("hi_zoomout()")
        writer.lua_model.append("hideconsole()")
        writer.save_as("heatflow_test.feh")
        writer.analyze()
        writer.load_solution()

        # select the material
        writer.lua_model.append("ho_selectblock(0, 0.1)")
        writer.lua_model.append("Fx, Fy = ho_blockintegral(3)")
        writer.write_out_result("Fx", "Fx")
        writer.write_out_result("Fy", "Fy")
        writer.close()

        try:
            reference = files("tests.integration_tests").joinpath("heatflow_test.lua")
            with open(reference) as f:
                content = f.readlines()
                counter_test = Counter(content)
                counter_reference = Counter(writer.lua_model)

                for key in counter_reference.keys():
                    # print(f'|{key}|', counter_reference[key.rstrip()], counter_test[key + "\n"])

                    # filter out path related commands
                    if "remove(" in key:
                        continue

                    if "saveas" in key:
                        continue

                    if "openfile" in key:
                        continue

                    self.assertEqual(
                        counter_reference[key.rstrip()],
                        counter_test[key + "\n"],
                    )

        except FileNotFoundError:
            self.assertTrue(False)
