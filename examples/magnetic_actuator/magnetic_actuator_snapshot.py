from adze_modeler.femm_wrapper import FemmExecutor, FemmWriter, femm_magnetic, MagneticMaterial, MagneticDirichlet
from adze_modeler.geometry import Geometry

"""
The description of the problem can be found here:
http://wiki.maxwell.sze.hu/index.php/Feladat_2

This script solves the 3. variant:
    - Excitation current: 2 A
    - Number of turns in the coil: 450
    - Used magnet: AlNiCo5 (mu_r = 19.97)

The force exerted on the magnet should be around -4.85 N.
"""

aluminium = MagneticMaterial("aluminium", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
air = MagneticMaterial("air", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
steel = MagneticMaterial("steel", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
magnet = MagneticMaterial("Alnico5", 19.97, 19.97, 51e3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


def create_snapshot():
    geo = Geometry()
    geo.import_svg("actuator.svg")

    # aligning the geometry to the axes
    minx = min(geo.nodes, key=lambda ni: ni.x).x
    miny = min(geo.nodes, key=lambda ni: ni.y).y
    maxx = max(geo.nodes, key=lambda ni: ni.x).x
    maxy = max(geo.nodes, key=lambda ni: ni.y).y
    offsetx = -minx  # align to the right side of y axis
    offsety = -miny - (maxy - miny) / 2  # align symmetrically to the x axis

    for i in range(len(geo.lines)):
        geo.lines[i].start_pt.x += offsetx
        geo.lines[i].start_pt.y += offsety

        geo.lines[i].end_pt.x += offsetx
        geo.lines[i].end_pt.y += offsety

    geo.generate_intersections()

    # create and save a model geometry in lua language
    femm_model = FemmWriter()
    femm_model.field = femm_magnetic
    femm_model.lua_model = femm_model.init_problem()
    femm_model.magnetic_problem(0, "millimeters", "axi", 1e-8)
    femm_model.lua_model += femm_model.create_geometry(geo)
    femm_model.lua_model.append("mi_zoomnatural()")
    femm_model.lua_model.append("hideconsole()")

    # Materials
    femm_model.add_material(air)
    femm_model.add_material(aluminium)
    femm_model.add_material(magnet)
    femm_model.add_material(steel)
    femm_model.add_bhcurve(steel, "bh_curve.dat", "h", "\t")

    # adding block labels
    femm_model.add_blocklabel(40, 0)  # Air
    femm_model.select_label(40, 0)
    femm_model.set_blockprop("air")
    femm_model.clear_selected()

    femm_model.add_blocklabel(5, -14)  # Lower aluminium
    femm_model.add_blocklabel(5, 14)  # Upper aluminium
    femm_model.select_label(5, -14)
    femm_model.select_label(5, 14)
    femm_model.set_blockprop("aluminium")
    femm_model.clear_selected()

    femm_model.add_blocklabel(3, 0)  # Magnet
    femm_model.select_label(3, 0)
    femm_model.set_blockprop("Alnico5", magdirection=90)
    femm_model.clear_selected()

    femm_model.add_blocklabel(14, 0)  # Steel
    femm_model.select_label(14, 0)
    femm_model.set_blockprop("steel")
    femm_model.clear_selected()

    I = 2.0  # Excitation current
    N = 450  # number of turns
    femm_model.add_circprop("Coil", I, 1)
    femm_model.add_blocklabel(11, 8)  # Upper Coil
    femm_model.select_label(11, 8)
    femm_model.set_blockprop("air", circuit_name="Coil", turns=-N)
    femm_model.clear_selected()

    femm_model.add_blocklabel(11, -8)  # Lower Coil
    femm_model.select_label(11, -8)
    femm_model.set_blockprop("air", circuit_name="Coil", turns=N)
    femm_model.clear_selected()

    # boundary conditions
    # dirichlet
    a0 = MagneticDirichlet("a0", 0, 0, 0, 0)  # "name", "a_0", "a_1", "a_2", "phi"
    femm_model.add_boundary(a0)

    # set dirichlet boundary condition on the bounding box
    femm_model.select_segment(35, 50)  # upper segment
    femm_model.select_segment(50, 0)  # right segment
    femm_model.select_segment(35, -50)  # lower segment

    femm_model.set_segment_prop("a0")
    femm_model.clear_selected()
    femm_model.lua_model.append("mi_refreshview()")

    femm_model.save_as("actuator.fem")
    femm_model.analyze()
    femm_model.load_solution()

    # Postprocessing
    # Calculate the force acting on the magnet
    femm_model.lua_model.append("mo_selectblock(0.1, 0)")
    femm_model.lua_model.append("Fy=mo_blockintegral(19)")  # manual: page 94
    femm_model.write_out_result("Fy", "Fy")  # write the results into femm_data.csv

    femm_model.close()  # close femm
    femm_model.write("actuator.lua")  # dump the entire script into actuator.lua


if __name__ == "__main__":
    create_snapshot()  # create the model
    FemmExecutor().run_femm("actuator.lua")  # run the script

    # Examine the results
    with open("femm_data.csv") as f:
        content = f.readline()
        Fy = float(content.strip().split(",")[1])

        print(f"Fy: {Fy:.4f} N")
