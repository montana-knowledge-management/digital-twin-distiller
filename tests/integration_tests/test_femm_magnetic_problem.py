import os
import unittest
from collections import Counter
from math import pi
from os import remove
from importlib_resources import files

from digital_twin_distiller.femm_wrapper import FemmWriter, MagneticMaterial, MagneticMixed


class TestFemmWriterWithExecutor(unittest.TestCase):
    """Tries to run a simple FEMM model, which is made with the FemmWriter class"""

    def test_air_core_coil_inductance_lua(self):
        """
        mi_loadsolution;
        c = mo_getcircuitproperties('icoil');
        y = c(3);
        closefemm;
        """
        # n = 100, ri = 1, ro = 2, z = 1
        # FEMM inductance =  0.64707 mH

        writer = FemmWriter()
        writer.lua_model = []  #
        writer.init_problem(out_file="femm_data.csv")

        # problem definition
        writer.magnetic_problem(0, "inches", "axi")

        # model geometry
        # rectangle (coil) -- ri, -z / 2, ro, z / 2 --
        n = 100
        ri = 1.0
        ro = 2.0
        z = 1.0
        r = 2.0 * max(ri, ro, z)

        # nodes
        z = 0.5 * z
        writer.add_node(ri, -z)
        writer.add_node(ri, z)
        writer.add_node(ro, -z)
        writer.add_node(ro, z)

        writer.add_segment(ri, -z, ri, z)
        writer.add_segment(ri, z, ro, z)
        writer.add_segment(ro, z, ro, -z)
        writer.add_segment(ro, -z, ri, -z)

        # air region
        writer.add_node(0, -r)
        writer.add_node(0, r)
        writer.add_segment(0, -r, 0, r)

        writer.add_arc(0.0, -r, 0.0, r, 180, 5)

        # circle properties and material definitions
        writer.add_circprop("icoil", 1, 1)

        writer.add_blocklabel((ri + ro) / 2, 0)
        writer.add_blocklabel(0.75 * r, 0)

        # add material properties
        coil = MagneticMaterial("coil", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
        air = MagneticMaterial("air", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)

        writer.add_material(coil)
        writer.add_material(air)

        # add boundary properties
        mixed_boundary = MagneticMixed("abc", 1.0 / (r * 0.0254 * pi * 4e-7), 0)
        writer.add_boundary(mixed_boundary)

        # set coil property
        writer.select_label((ri + ro) / 2, 0)
        writer.set_blockprop("coil", 0, r / 20, 0, circuit_name="icoil", turns=n, magdirection=0)
        writer.clear_selected()

        # set air
        writer.select_label(0.75 * r, 0)
        writer.set_blockprop("air", 0, r / 100, 0, circuit_name="<None>", turns=0, magdirection=0)
        writer.clear_selected()

        # set boundaries
        writer.select_arc_segment(r, 0)
        writer.set_arc_segment_prop(5, "abc", 0, 0)

        # the model have to be saved into a femm format to be used from the FEMM-wrapper script
        writer.save_as("test.fem")
        writer.analyze()
        writer.load_solution()
        writer.get_circuit_properties("icoil", result="current, volt, flux")
        writer.write_out_result("current", "current")
        writer.write_out_result("volt", "volt")
        writer.write_out_result("flux", "flux")

        writer.close()
        writer.write("magnetic_ref.lua")

        try:
            reference = files("tests.integration_tests").joinpath("magnetic.lua")
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

            os.remove('magnetic_ref.lua')

        except FileNotFoundError:
            self.assertTrue(False)
