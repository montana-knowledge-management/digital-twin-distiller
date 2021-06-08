import os

from adze_modeler.femm_wrapper import FemmExecutor
from adze_modeler.femm_wrapper import FemmWriter
from adze_modeler.femm_wrapper import HeatFlowConvection
from adze_modeler.femm_wrapper import HeatFlowFixedTemperature
from adze_modeler.femm_wrapper import HeatFlowMaterial
from adze_modeler.femm_wrapper import femm_heat_flow


def c2k(C):
    return C + 274.15


def test_heat_problem():
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
    writer.write("heatflow_test.lua")


def run():
    test_heat_problem()
    FemmExecutor().run_femm("heatflow_test.lua")

    with open("heat_data.csv") as f:
        content = f.readlines()
        # print(content[0])
        # print(content[1])
        Fx = float(content[0].split(",")[1])
        Fy = float(content[1].split(",")[1])

        os.remove("heat_data.csv")
        os.remove("heatflow_test.anh")
        os.remove("heatflow_test.feh")
        os.remove("heatflow_test.lua")

        # expected results
        exp_Fx = 0.0112
        exp_Fy = 9551.0549

        if (exp_Fx - round(Fx, 4) + exp_Fy - round(Fy, 4)) < 1e-3:
            return True

        return False


if __name__ == "__main__":
    run()
