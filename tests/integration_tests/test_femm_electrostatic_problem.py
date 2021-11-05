import unittest
from collections import Counter

from importlib_resources import files

from adze_modeler.femm_wrapper import ElectrostaticMaterial, FemmWriter, femm_electrostatic


class TestFemmElectrostaticProblem(unittest.TestCase):
    # integration test ignored from the unittest list
    def test_electrostatic_problem(self):
        writer = FemmWriter()
        writer.field = femm_electrostatic
        writer.init_problem("electrostatic_data.csv")

        writer.electrostatic_problem("centimeters", "planar")

        a = 10  # cm
        epsilon_r = 2.1  # -
        Ug = 10  # V

        writer.add_node(-a / 2, a / 2)
        writer.add_node(-a / 2, 0)
        writer.add_node(0, 0)
        writer.add_node(0, -a / 2)
        writer.add_node(a / 2, -a / 2)
        writer.add_node(a / 2, a / 2)

        writer.add_segment(-a / 2, a / 2, -a / 2, 0)
        writer.add_segment(-a / 2, 0, 0, 0)
        writer.add_segment(0, 0, 0, -a / 2)
        writer.add_segment(0, -a / 2, a / 2, -a / 2)
        writer.add_segment(a / 2, -a / 2, a / 2, a / 2)
        writer.add_segment(a / 2, a / 2, -a / 2, a / 2)

        # Adding material properties
        blocklabel = (a / 4, a / 4)
        mat = ElectrostaticMaterial("Teflon", epsilon_r, epsilon_r, 0)
        writer.add_material(mat)
        writer.add_blocklabel(*blocklabel)
        writer.select_label(*blocklabel)
        writer.set_blockprop("Teflon")

        # Adding boundary properties
        writer.add_pointprop("Ug", Vp=Ug)
        writer.add_pointprop("U0", Vp=0)

        writer.select_node(0, 0)
        writer.set_pointprop("Ug")
        writer.lua_model.append("ei_clearselected()")

        writer.select_node(a / 2, a / 2)
        writer.set_pointprop("U0")
        writer.lua_model.append("ei_clearselected()")

        writer.lua_model.append("ei_zoomnatural()")
        writer.lua_model.append("ei_zoomout()")
        writer.lua_model.append("hideconsole()")
        writer.save_as("electrostatic_test.fee")
        writer.analyze()
        writer.load_solution()

        # Examine the results
        writer.lua_model.append(f"eo_selectblock({a / 4}, {a / 4})")
        writer.lua_model.append("E = eo_blockintegral(0)")  # Stored Energy
        writer.write_out_result("E", "E")
        writer.close()

        try:
            reference = files("tests.integration_tests").joinpath("electrostatic_test.lua")
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
