import os

from adze_modeler.femm_wrapper import CurrentFlowFixedVoltage
from adze_modeler.femm_wrapper import CurrentFlowMaterial
from adze_modeler.femm_wrapper import CurrentFlowSurfaceCurrent
from adze_modeler.femm_wrapper import FemmExecutor
from adze_modeler.femm_wrapper import FemmWriter
from adze_modeler.femm_wrapper import kw_current_flow

# Problem parameters

wt = 85
wk = 56.6
wcm = 38
ht = 36
hk = 7.6
hcm = 13
hc = 10
dk = 8.2

Ig = 400
depth = 5
J = Ig / (depth * 1e-3 * ht * 1e-3)


def current_flow_problem():
    writer = FemmWriter()
    writer.field = kw_current_flow
    writer.init_problem("current_data.csv")

    writer.currentflow_problem("millimeters", "planar", depth=depth)

    # Geometry points
    writer.add_node(-wcm / 2, hc / 2)
    writer.add_node(wcm / 2, hc / 2)
    writer.add_node(-wcm / 2, -hc / 2)
    writer.add_node(wcm / 2, -hc / 2)

    writer.add_node(-wcm / 2, hc / 2 + hcm)
    writer.add_node(wcm / 2, hc / 2 + hcm)
    writer.add_node(-wcm / 2, -hc / 2 - hcm)
    writer.add_node(wcm / 2, -hc / 2 - hcm)

    writer.add_node(-wt / 2, ht / 2)
    writer.add_node(wt / 2, ht / 2)
    writer.add_node(-wt / 2, -ht / 2)
    writer.add_node(wt / 2, -ht / 2)

    # Vertical line segments
    writer.add_segment(-wt / 2, ht / 2, -wt / 2, -ht / 2)
    writer.add_segment(-wcm / 2, ht / 2, -wcm / 2, hc / 2)
    writer.add_segment(-wcm / 2, ht / 2, -wcm / 2, -hc / 2)
    writer.add_segment(-wcm / 2, -hc / 2, -wcm / 2, -ht / 2)

    writer.add_segment(wt / 2, ht / 2, wt / 2, -ht / 2)
    writer.add_segment(wcm / 2, ht / 2, wcm / 2, hc / 2)
    writer.add_segment(wcm / 2, ht / 2, wcm / 2, -hc / 2)
    writer.add_segment(wcm / 2, -hc / 2, wcm / 2, -ht / 2)

    # Horizontal line segments
    writer.add_segment(-wt / 2, ht / 2, -wcm / 2, ht / 2)
    writer.add_segment(-wcm / 2, ht / 2, wcm / 2, ht / 2)
    writer.add_segment(wcm / 2, ht / 2, wt / 2, ht / 2)

    writer.add_segment(-wt / 2, -ht / 2, -wcm / 2, -ht / 2)
    writer.add_segment(-wcm / 2, -ht / 2, wcm / 2, -ht / 2)
    writer.add_segment(wcm / 2, -ht / 2, wt / 2, -ht / 2)

    writer.add_segment(-wcm / 2, hc / 2, wcm / 2, hc / 2)
    writer.add_segment(-wcm / 2, -hc / 2, wcm / 2, -hc / 2)

    # Circle points
    writer.add_node(-wk / 2, hk / 2 + dk / 2)
    writer.add_node(-wk / 2 - dk, hk / 2 + dk / 2)

    writer.add_node(-wk / 2, -hk / 2 - dk / 2)
    writer.add_node(-wk / 2 - dk, -hk / 2 - dk / 2)

    writer.add_node(wk / 2, hk / 2 + dk / 2)
    writer.add_node(wk / 2 + dk, hk / 2 + dk / 2)

    writer.add_node(wk / 2, -hk / 2 - dk / 2)
    writer.add_node(wk / 2 + dk, -hk / 2 - dk / 2)

    # Adding arc segments
    writer.add_arc(-wk / 2, hk / 2 + dk / 2, -wk / 2 - dk, hk / 2 + dk / 2, 180, 1)
    writer.add_arc(-wk / 2 - dk, hk / 2 + dk / 2, -wk / 2, hk / 2 + dk / 2, 180, 1)

    writer.add_arc(-wk / 2, -hk / 2 - dk / 2, -wk / 2 - dk, -hk / 2 - dk / 2, 180, 1)
    writer.add_arc(-wk / 2 - dk, -hk / 2 - dk / 2, -wk / 2, -hk / 2 - dk / 2, 180, 1)

    writer.add_arc(wk / 2, hk / 2 + dk / 2, wk / 2 + dk, hk / 2 + dk / 2, 180, 1)
    writer.add_arc(wk / 2 + dk, hk / 2 + dk / 2, wk / 2, hk / 2 + dk / 2, 180, 1)

    writer.add_arc(wk / 2, -hk / 2 - dk / 2, wk / 2 + dk, -hk / 2 - dk / 2, 180, 1)
    writer.add_arc(wk / 2 + dk, -hk / 2 - dk / 2, wk / 2, -hk / 2 - dk / 2, 180, 1)

    # Adding block labels
    bl_copper = ((wt - wcm) / 4 + wcm / 2, 0)
    bl_cmanganin = (0, hc / 2 + hcm / 2)
    bl_titanium = (wk / 2 + dk / 2, hk / 2 + dk / 2)
    bl_nomesh = (0, 0)

    writer.add_blocklabel(bl_copper[0], bl_copper[1])
    writer.add_blocklabel(-bl_copper[0], bl_copper[1])

    writer.add_blocklabel(bl_cmanganin[0], bl_cmanganin[1])
    writer.add_blocklabel(bl_cmanganin[0], -bl_cmanganin[1])

    writer.add_blocklabel(bl_titanium[0], bl_titanium[1])
    writer.add_blocklabel(bl_titanium[0], -bl_titanium[1])
    writer.add_blocklabel(-bl_titanium[0], bl_titanium[1])
    writer.add_blocklabel(-bl_titanium[0], -bl_titanium[1])

    writer.add_blocklabel(bl_nomesh[0], bl_nomesh[1])

    # Adding materials
    mat_copper = CurrentFlowMaterial("Copper", 58e6, 58e6, 0, 0, 0, 0)
    mat_cmanganin = CurrentFlowMaterial("Copper-Manganin", 20.833e6, 20.833e6, 0, 0, 0, 0)
    mat_titanium = CurrentFlowMaterial("Titanium", 1.789e6, 1.789e6, 0, 0, 0, 0)

    writer.add_material(mat_copper)
    writer.add_material(mat_cmanganin)
    writer.add_material(mat_titanium)

    writer.select_label(bl_copper[0], bl_copper[1])
    writer.select_label(-bl_copper[0], bl_copper[1])
    writer.set_blockprop("Copper")
    writer.lua_model.append("ci_clearselected()")

    writer.select_label(bl_cmanganin[0], bl_cmanganin[1])
    writer.select_label(bl_cmanganin[0], -bl_cmanganin[1])
    writer.set_blockprop("Copper-Manganin")
    writer.lua_model.append("ci_clearselected()")

    writer.select_label(bl_titanium[0], bl_titanium[1])
    writer.select_label(bl_titanium[0], -bl_titanium[1])
    writer.select_label(-bl_titanium[0], bl_titanium[1])
    writer.select_label(-bl_titanium[0], -bl_titanium[1])
    writer.set_blockprop("Titanium")
    writer.lua_model.append("ci_clearselected()")

    writer.select_label(bl_nomesh[0], bl_nomesh[1])
    writer.set_blockprop("<No Mesh>")
    writer.lua_model.append("ci_clearselected()")

    # Adding boundary properties
    excitation = CurrentFlowSurfaceCurrent("Jin", J)
    ground = CurrentFlowFixedVoltage("GND", 0)

    writer.add_boundary(excitation)
    writer.add_boundary(ground)

    writer.select_segment(-wt / 2, 0)
    writer.set_segment_prop("Jin")
    writer.lua_model.append("ci_clearselected()")

    writer.select_segment(wt / 2, 0)
    writer.set_segment_prop("GND")
    writer.lua_model.append("ci_clearselected()")

    writer.lua_model.append("ci_zoomnatural()")
    writer.lua_model.append("ci_zoomout()")
    writer.lua_model.append("hideconsole()")
    writer.save_as("current_test.fec")
    writer.analyze()
    writer.load_solution()

    # Examine the results
    writer.lua_model.append(f"co_selectblock({bl_copper[0]}, {bl_copper[1]})")
    writer.lua_model.append(f"co_selectblock({-bl_copper[0]}, {bl_copper[1]})")

    writer.lua_model.append(f"co_selectblock({bl_cmanganin[0]}, {bl_cmanganin[1]})")
    writer.lua_model.append(f"co_selectblock({bl_cmanganin[0]}, {-bl_cmanganin[1]})")

    writer.lua_model.append(f"co_selectblock({bl_titanium[0]}, {bl_titanium[1]})")
    writer.lua_model.append(f"co_selectblock({bl_titanium[0]}, {-bl_titanium[1]})")
    writer.lua_model.append(f"co_selectblock({-bl_titanium[0]}, {bl_titanium[1]})")
    writer.lua_model.append(f"co_selectblock({-bl_titanium[0]}, {-bl_titanium[1]})")

    writer.lua_model.append("P = co_blockintegral(0)")  # Power Loss
    writer.lua_model.append(writer.write_out_result("P", "P"))
    writer.close()

    # Execute the script
    writer.write("current_test.lua")


def run():
    current_flow_problem()
    FemmExecutor().run_femm("current_test.lua")

    with open("current_data.csv") as f:
        content = f.readlines()
        P = float(content[0].split(",")[1])  # Power loss
        R = P / Ig ** 2  # Resistance should be around 20 uOhm

        expected_result = round(2.0069e-5, 4)
        result = round(float(R), 4)

        try:
            os.remove("current_data.csv")
            os.remove("current_test.anc")
            os.remove("current_test.fec")
            os.remove("current_test.lua")
        except FileNotFoundError:
            print("Files not generated")
            return False

        if expected_result == result:
            return True

    return False


if __name__ == "__main__":
    run()
