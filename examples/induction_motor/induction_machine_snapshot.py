from adze_modeler.geometry import Geometry
from adze_modeler.objects import Node
from adze_modeler.femm_wrapper import FemmWriter, MagneticMaterial, femm_magnetic, MagneticDirichlet, MagneticAnti
from adze_modeler.material import MaterialSnapshot
from adze_modeler.boundaries import BoundarySnaphot
from adze_modeler.femm_wrapper import FemmExecutor
from math import pi

# material properties
# silicon core with linea bh relationships
silicon_core = MagneticMaterial("silicon_core", 7000, 7000, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
aluminum = MagneticMaterial("aluminum", 1, 1, 0, 0, 34.5, 0, 0, 1, 0, 0, 0, 0, 0)
coil = MagneticMaterial("copper", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
air = MagneticMaterial("air", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)


def create_phase(base_node, group, circuit_name, phase_current):
    winding = MaterialSnapshot()
    winding.material_definition = coil
    winding.field_type = femm_magnetic  # magnetic field with the femm material

    winding.turn_umber = 44
    winding.group = group
    winding.circuit_name = circuit_name
    winding.circuit_current = phase_current
    winding.circuit_type = 1  # series connected windings

    cmds = []
    for i in range(0, 3):
        winding.region_labels.append(base_node)
        base_node = base_node.rotate(rad=10.0 / 180 * pi)
        cmds = winding.create_femm_block_label()

    return cmds


def create_snapshot():
    """Create a snapshot manually, where every quantity can be calculated for one given position of the motor."""

    # the geometry is defined in a dxf file
    geo = Geometry()
    geo.import_dxf("2horse.dxf")

    # create and save a model geometry in lua language
    femm_model = FemmWriter()
    femm_model.lua_model = femm_model.init_problem()
    femm_model.lua_model.append(FemmWriter().magnetic_problem(4, "millimeters", "planar", 1e-8, 100, 20))
    femm_model.lua_model += femm_model.create_geometry(geo)

    # airgap ----------------------------------------------------------------------
    airgap = MaterialSnapshot()
    airgap.material_definition = air
    airgap.field_type = femm_magnetic  # magnetic field with the femm material

    # label positions
    airgap.region_labels = [Node(40.5, 3.5)]
    femm_model.lua_model += airgap.create_femm_block_label()

    # cores -----------------------------------------------------------------------
    cores = MaterialSnapshot()
    cores.material_definition = silicon_core
    cores.field_type = femm_magnetic  # magnetic field with the femm material

    # label positions
    cores.region_labels = [Node(16.0, 16.0), Node(40.0, 40.0)]
    femm_model.lua_model += cores.create_femm_block_label()

    # aluminum rods ----------------------------------------------------------------
    rods = MaterialSnapshot()
    rods.material_definition = aluminum
    rods.field_type = femm_magnetic  # magnetic field with the femm material

    base_node = Node(31.6, 0.33)
    for i in range(0, 8):
        rods.region_labels.append(base_node)
        base_node = base_node.rotate(rad=12.5 / 180 * pi)
    femm_model.lua_model += rods.create_femm_block_label()

    # coil definitions
    femm_model.lua_model += create_phase(Node(45.5, 3.5), 1, 'A', '1')
    femm_model.lua_model += create_phase(Node(20.0, 43.9), 3, 'B', '-0.5+I*0.8660254037844386')
    femm_model.lua_model += create_phase(Node(39.285, 26.8865), 2, 'C', '-0.5-I*0.866025403784439')

    # boundary conditions
    # dirichlet
    a0 = MagneticDirichlet('a0', 0, 0, 0, 0)  # "name", "a_0", "a_1", "a_2", "phi"

    diri = BoundarySnaphot()
    diri.boundary = a0
    diri.group = 1
    diri.field_type = femm_magnetic  # magnetic field with the femm material

    for arc in geo.circle_arcs:
        if abs((float(arc.start_pt.x) ** 2 + float(arc.start_pt.y) ** 2) ** 0.5 - 12.5) < 1e-3:
            diri.elements.append(arc)
        if abs((float(arc.start_pt.x) ** 2 + float(arc.start_pt.y) ** 2) ** 0.5 - 65.0) < 1e-3:
            diri.elements.append(arc)
    # diri.elements = []
    femm_model.lua_model += diri.create_femm_boundaries()

    # antiperiodic conditions
    ap1 = MagneticAnti('ap1')
    app1 = BoundarySnaphot()
    app1.boundary = ap1
    app1.group = 2
    app1.field_type = femm_magnetic  # magnetic field with the femm material

    ap2 = MagneticAnti('ap2')
    app2 = BoundarySnaphot()
    app2.boundary = ap2
    app2.group = 3
    app2.field_type = femm_magnetic  # magnetic field with the femm material

    ap3 = MagneticAnti('ap3')
    app3 = BoundarySnaphot()
    app3.boundary = ap3
    app3.group = 4
    app3.field_type = femm_magnetic  # magnetic field with the femm material

    ap4 = MagneticAnti('ap4')
    app4 = BoundarySnaphot()
    app4.boundary = ap4
    app4.group = 5
    app4.field_type = femm_magnetic  # magnetic field with the femm material

    for line in geo.lines:
        ## 1
        if line.start_pt.y < 1e-3:
            if line.start_pt.x < 25.001 and line.end_pt.x < 25.001:
                app1.elements.append(line)

        if line.start_pt.x < 1e-3:
            if line.start_pt.y < 25.001 and line.end_pt.y < 25.001:
                app1.elements.append(line)

        ## 2
        if line.start_pt.y < 1e-3:
            if line.start_pt.x < 39.501 and line.end_pt.x < 39.501 and line.start_pt.x > 25.0:
                app2.elements.append(line)

        if line.start_pt.x < 1e-3:
            if line.start_pt.y < 39.501 and line.end_pt.y < 39.501 and line.start_pt.y > 25.0:
                app2.elements.append(line)

        ## 3
        if line.start_pt.y < 1e-3:
            if line.start_pt.x < 40.376 and line.end_pt.x < 40.376 and line.start_pt.x > 39.50:
                app3.elements.append(line)

        if line.start_pt.x < 1e-3:
            if line.start_pt.y < 40.376 and line.end_pt.y < 40.376 and line.start_pt.y > 39.50:
                app3.elements.append(line)

        ## 4
        if line.start_pt.y < 1e-3:
            if line.start_pt.x < 65.1 and line.end_pt.x < 65.1 and line.start_pt.x > 40.3 and line.end_pt.x > 40.3:
                app4.elements.append(line)

        if line.start_pt.x < 1e-3:
            if line.start_pt.y < 65.1 and line.end_pt.y < 65.1 and line.start_pt.y > 40.3 and line.end_pt.y > 40.3:
                app4.elements.append(line)

    femm_model.lua_model += app1.create_femm_boundaries()
    femm_model.lua_model += app2.create_femm_boundaries()
    femm_model.lua_model += app3.create_femm_boundaries()
    femm_model.lua_model += app4.create_femm_boundaries()

    # save to file
    femm_model.lua_model.extend(femm_model.close())
    femm_model.write('2horse.lua')
    return


if __name__ == "__main__":
    create_snapshot()
    FemmExecutor().run_femm('2horse.lua')

