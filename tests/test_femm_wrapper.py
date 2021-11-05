import os
import warnings
from pathlib import Path
from unittest import TestCase

import adze_modeler.objects as obj
from adze_modeler.femm_wrapper import (
    CurrentFlowAntiPeriodic,
    CurrentFlowFixedVoltage,
    CurrentFlowMaterial,
    CurrentFlowMixed,
    CurrentFlowPeriodic,
    CurrentFlowSurfaceCurrent,
    ElectrostaticAntiPeriodic,
    ElectrostaticFixedVoltage,
    ElectrostaticMaterial,
    ElectrostaticMixed,
    ElectrostaticPeriodic,
    ElectrostaticSurfaceCharge,
    FemmExecutor,
    FemmWriter,
    HeatFlowAntiPeriodic,
    HeatFlowConvection,
    HeatFlowFixedTemperature,
    HeatFlowHeatFlux,
    HeatFlowMaterial,
    HeatFlowPeriodic,
    HeatFlowRadiation,
    MagneticDirichlet,
    MagneticMaterial,
    MagneticMixed,
    femm_current_flow,
    femm_electrostatic,
    femm_heat_flow,
    femm_magnetic,
)
from adze_modeler.geometry import Geometry


class FemmTester(TestCase):
    def test_not_defined_writer(self):
        writer = FemmWriter()
        writer.field = None

        self.assertRaises(ValueError)

    def test_validate_field(self):
        writer = FemmWriter()
        writer.field = femm_electrostatic
        self.assertEqual(writer.validate_field(), True)

        writer.field = femm_heat_flow
        self.assertEqual(writer.validate_field(), True)

        writer.field = femm_current_flow
        self.assertEqual(writer.validate_field(), True)

        writer.field = femm_magnetic
        self.assertEqual(writer.validate_field(), True)

        writer.field = "alma"
        self.assertRaises(ValueError, writer.validate_field)

        writer.field = "eper"
        self.assertRaises(ValueError, writer.validate_field, "eper")

        writer.field = femm_electrostatic
        self.assertRaises(ValueError, writer.validate_field, "eper")

    def test_validate_units(self):
        writer = FemmWriter()
        self.assertEqual(True, writer.validate_units("inches"))
        self.assertEqual(True, writer.validate_units("millimeters"))
        self.assertEqual(True, writer.validate_units("centimeters"))
        self.assertEqual(True, writer.validate_units("mils"))
        self.assertEqual(True, writer.validate_units("meters"))
        self.assertEqual(True, writer.validate_units("micrometers"))
        self.assertRaises(ValueError, writer.validate_units, "alma")

    def test_write(self):
        writer = FemmWriter()
        writer.field = femm_heat_flow
        writer.lua_model.append("alma")
        writer.write("test_write.lua")
        self.assertEqual(True, os.path.exists("test_write.lua"))
        os.remove("test_write.lua")

    def test_addnode(self):
        x = 1.0
        y = 0.0

        # magnetic field
        res = FemmWriter().add_node(x, y)
        self.assertEqual("mi_addnode(1.0, 0.0)", res)

        # current flow
        fmw = FemmWriter()
        fmw.field = femm_current_flow
        res = fmw.add_node(x, y)
        self.assertEqual("ci_addnode(1.0, 0.0)", res)

        # electrostatic
        fmw = FemmWriter()
        fmw.field = femm_electrostatic
        res = fmw.add_node(x, y)
        self.assertEqual("ei_addnode(1.0, 0.0)", res)

        fmw = FemmWriter()
        fmw.field = femm_heat_flow
        res = fmw.add_node(x, y)
        self.assertEqual("hi_addnode(1.0, 0.0)", res)

    def test_add_segment(self):
        x1 = 1.0
        y1 = 0.0

        x2 = 1.0
        y2 = 1.0

        res = FemmWriter().add_segment(x1, y1, x2, y2)
        self.assertEqual("mi_addsegment(1.0, 0.0, 1.0, 1.0)", res)

        # current field
        fmw = FemmWriter()
        fmw.field = femm_current_flow
        res = fmw.add_segment(x1, y1, x2, y2)
        self.assertEqual("ci_addsegment(1.0, 0.0, 1.0, 1.0)", res)

        # electrostatic field
        fmw = FemmWriter()
        fmw.field = femm_electrostatic
        res = fmw.add_segment(x1, y1, x2, y2)
        self.assertEqual("ei_addsegment(1.0, 0.0, 1.0, 1.0)", res)

        # heat flow
        fmw = FemmWriter()
        fmw.field = femm_heat_flow
        res = fmw.add_segment(x1, y1, x2, y2)
        self.assertEqual("hi_addsegment(1.0, 0.0, 1.0, 1.0)", res)

    def test_addblocklabel(self):
        x = 1.0
        y = 0.0

        res = FemmWriter().add_blocklabel(x, y)
        self.assertEqual("mi_addblocklabel(1.0, 0.0)", res)

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        res = fmw.add_blocklabel(x, y)
        self.assertEqual("ci_addblocklabel(1.0, 0.0)", res)

        fmw.field = femm_heat_flow
        res = fmw.add_blocklabel(x, y)
        self.assertEqual("hi_addblocklabel(1.0, 0.0)", res)

        fmw.field = femm_electrostatic
        res = fmw.add_blocklabel(x, y)
        self.assertEqual("ei_addblocklabel(1.0, 0.0)", res)

    def test_addarc(self):
        x1, x2 = 1.0, 1.0
        y1, y2 = 0.0, 1.0

        res = FemmWriter().add_arc(x1, y1, x2, y2, 90.0, 1)
        self.assertEqual("mi_addarc(1.0, 0.0, 1.0, 1.0, 90.0, 1)", res)

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        res = fmw.add_arc(x1, y1, x2, y2, 90.0, 1)
        self.assertEqual("ci_addarc(1.0, 0.0, 1.0, 1.0, 90.0, 1)", res)

        fmw.field = femm_electrostatic
        res = fmw.add_arc(x1, y1, x2, y2, 90.0, 1)
        self.assertEqual("ei_addarc(1.0, 0.0, 1.0, 1.0, 90.0, 1)", res)

        fmw.field = femm_heat_flow
        res = fmw.add_arc(x1, y1, x2, y2, 90.0, 1)
        self.assertEqual("hi_addarc(1.0, 0.0, 1.0, 1.0, 90.0, 1)", res)

    def test_delete_selected(self):
        self.assertEqual("mi_deleteselected", FemmWriter().delete_selected())

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_deleteselected", fmw.delete_selected())

        fmw.field = femm_heat_flow
        self.assertEqual("hi_deleteselected", fmw.delete_selected())

        fmw.field = femm_electrostatic
        self.assertEqual("ei_deleteselected", fmw.delete_selected())

    def test_delete_selected_nodes(self):
        self.assertEqual("mi_deleteselectednodes", FemmWriter().delete_selected_nodes())

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_deleteselectednodes", fmw.delete_selected_nodes())

        fmw.field = femm_heat_flow
        self.assertEqual("hi_deleteselectednodes", fmw.delete_selected_nodes())

        fmw.field = femm_electrostatic
        self.assertEqual("ei_deleteselectednodes", fmw.delete_selected_nodes())

    def test_delete_selected_labels(self):
        self.assertEqual("mi_deleteselectedlabels", FemmWriter().delete_selected_labels())

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_deleteselectedlabels", fmw.delete_selected_labels())

        fmw.field = femm_heat_flow
        self.assertEqual("hi_deleteselectedlabels", fmw.delete_selected_labels())

        fmw.field = femm_electrostatic
        self.assertEqual("ei_deleteselectedlabels", fmw.delete_selected_labels())

    def test_delete_selected_segments(self):
        self.assertEqual("mi_deleteselectedsegments", FemmWriter().delete_selected_segments())

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_deleteselectedsegments", fmw.delete_selected_segments())

        fmw.field = femm_heat_flow
        self.assertEqual("hi_deleteselectedsegments", fmw.delete_selected_segments())

        fmw.field = femm_electrostatic
        self.assertEqual("ei_deleteselectedsegments", fmw.delete_selected_segments())

    def test_delete_selected_arc_segments(self):
        self.assertEqual(
            "mi_deleteselectedarcsegments",
            FemmWriter().delete_selected_arc_segments(),
        )

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_deleteselectedarcsegments", fmw.delete_selected_arc_segments())

        fmw.field = femm_heat_flow
        self.assertEqual("hi_deleteselectedarcsegments", fmw.delete_selected_arc_segments())

        fmw.field = femm_electrostatic
        self.assertEqual("ei_deleteselectedarcsegments", fmw.delete_selected_arc_segments())

    def test_clear_seelcted(self):
        self.assertEqual("mi_clearselected()", FemmWriter().clear_selected())

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_clearselected()", fmw.clear_selected())

        fmw.field = femm_heat_flow
        self.assertEqual("hi_clearselected()", fmw.clear_selected())

        fmw.field = femm_electrostatic
        self.assertEqual("ei_clearselected()", fmw.clear_selected())

    def test_select_segment(self):
        self.assertEqual("mi_selectsegment(1.0, 1.0)", FemmWriter().select_segment(1.0, 1.0))

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_selectsegment(1.0, 1.0)", fmw.select_segment(1.0, 1.0))

        fmw.field = femm_heat_flow
        self.assertEqual("hi_selectsegment(1.0, 1.0)", fmw.select_segment(1.0, 1.0))

        fmw.field = femm_electrostatic
        self.assertEqual("ei_selectsegment(1.0, 1.0)", fmw.select_segment(1.0, 1.0))

    def test_select_node(self):
        self.assertEqual("mi_selectnode(1.0, 1.0)", FemmWriter().select_node(1.0, 1.0))

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_selectnode(1.0, 1.0)", fmw.select_node(1.0, 1.0))

        fmw.field = femm_heat_flow
        self.assertEqual("hi_selectnode(1.0, 1.0)", fmw.select_node(1.0, 1.0))

        fmw.field = femm_electrostatic
        self.assertEqual("ei_selectnode(1.0, 1.0)", fmw.select_node(1.0, 1.0))

    def test_select_label(self):
        self.assertEqual("mi_selectlabel(1.0, 1.0)", FemmWriter().select_label(1.0, 1.0))

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_selectlabel(1.0, 1.0)", fmw.select_label(1.0, 1.0))

        fmw.field = femm_heat_flow
        self.assertEqual("hi_selectlabel(1.0, 1.0)", fmw.select_label(1.0, 1.0))

        fmw.field = femm_electrostatic
        self.assertEqual("ei_selectlabel(1.0, 1.0)", fmw.select_label(1.0, 1.0))

    def test_select_group(self):
        self.assertEqual("mi_selectgroup(4)", FemmWriter().select_group(4))

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual("ci_selectgroup(4)", fmw.select_group(4))

        fmw.field = femm_heat_flow
        self.assertEqual("hi_selectgroup(4)", fmw.select_group(4))

        fmw.field = femm_electrostatic
        self.assertEqual("ei_selectgroup(4)", fmw.select_group(4))

    def test_select_circle(self):
        self.assertEqual(
            "mi_selectcircle(1.0, 2.0, 0.4, 3)",
            FemmWriter().select_circle(1.0, 2.0, 0.4, 3),
        )

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual(
            "ci_selectcircle(1.0, 2.0, 0.4, 3)",
            fmw.select_circle(1.0, 2.0, 0.4, 3),
        )

        fmw.field = femm_heat_flow
        self.assertEqual(
            "hi_selectcircle(1.0, 2.0, 0.4, 3)",
            fmw.select_circle(1.0, 2.0, 0.4, 3),
        )

        fmw.field = femm_electrostatic
        self.assertEqual(
            "ei_selectcircle(1.0, 2.0, 0.4, 3)",
            fmw.select_circle(1.0, 2.0, 0.4, 3),
        )

    def test_select_rectangle(self):
        self.assertEqual(
            "mi_selectrectangle(1.0,2.0,3.0,4.0,3)",
            FemmWriter().select_rectangle(1.0, 2.0, 3.0, 4.0, 3),
        )

        fmw = FemmWriter()
        fmw.field = femm_current_flow
        self.assertEqual(
            "ci_selectrectangle(1.0,2.0,3.0,4.0,3)",
            fmw.select_rectangle(1.0, 2.0, 3.0, 4.0, 3),
        )

        fmw.field = femm_heat_flow
        self.assertEqual(
            "hi_selectrectangle(1.0,2.0,3.0,4.0,3)",
            fmw.select_rectangle(1.0, 2.0, 3.0, 4.0, 3),
        )

        fmw.field = femm_electrostatic
        self.assertEqual(
            "ei_selectrectangle(1.0,2.0,3.0,4.0,3)",
            fmw.select_rectangle(1.0, 2.0, 3.0, 4.0, 3),
        )

    def test_magnetic_problem(self):
        self.assertEqual(
            r"mi_probdef(50,'millimeters','axi',1e-08, 1, 30, 0)",
            FemmWriter().magnetic_problem(50, "millimeters", "axi"),
        )

    def test_heat_problem(self):
        writer = FemmWriter()
        writer.field = femm_heat_flow
        self.assertEqual(
            'hi_probdef("inches", "planar", 1e-08, 1, 30, "", 0)',
            writer.heat_problem("inches", "planar"),
        )
        self.assertRaises(ValueError, writer.heat_problem, "eper", "planar")
        self.assertRaises(ValueError, writer.heat_problem, "meters", "qwertz")

    def test_current_flow_problem(self):
        writer = FemmWriter()
        writer.field = femm_current_flow
        self.assertEqual(
            'ci_probdef("inches", "planar", 100, 1e-08, 1, 30)',
            writer.currentflow_problem("inches", "planar", 100, 1e-8, 1, 30),
        )
        self.assertRaises(ValueError, writer.currentflow_problem, "eper", "planar")
        self.assertRaises(ValueError, writer.currentflow_problem, "meters", "qwertz")

    def test_electrostatic_problem(self):
        writer = FemmWriter()
        writer.field = femm_electrostatic
        self.assertEqual(
            'ei_probdef("inches", "planar", 1e-08, 1, 30)',
            writer.electrostatic_problem("inches", "planar", 1e-8, 1, 30),
        )
        self.assertRaises(ValueError, writer.electrostatic_problem, "barack", "planar")
        self.assertRaises(ValueError, writer.electrostatic_problem, "mils", "planadawdawr")

    def test_init_problem(self):
        # TODO: check it
        # self.assertIn("showconsole", FemmWriter().init_problem()[1])
        # self.assertIn("clear", FemmWriter().init_problem()[1])

        writer = FemmWriter()
        writer.field = femm_magnetic
        print(writer.init_problem())
        self.assertEqual("newdocument(0)", writer.init_problem()[1])

        writer.field = femm_electrostatic
        self.assertEqual("newdocument(1)", writer.init_problem()[1])

        writer.field = femm_heat_flow
        self.assertEqual("newdocument(2)", writer.init_problem()[1])

        writer.field = femm_current_flow
        self.assertEqual("newdocument(3)", writer.init_problem()[1])

    def test_close(self):
        writer = FemmWriter()
        writer.field = femm_electrostatic
        self.assertEqual("eo_close()", writer.close()[1])
        self.assertEqual("ei_close()", writer.close()[2])

        writer.field = femm_heat_flow
        self.assertEqual("ho_close()", writer.close()[1])
        self.assertEqual("hi_close()", writer.close()[2])

        writer.field = femm_current_flow
        self.assertEqual("co_close()", writer.close()[1])
        self.assertEqual("ci_close()", writer.close()[2])

        writer.field = femm_magnetic
        self.assertEqual("mo_close()", writer.close()[1])
        self.assertEqual("mi_close()", writer.close()[2])

    def test_add_circ_prop(self):
        self.assertEqual(
            'mi_addcircprop("test",1,0)',
            FemmWriter().add_circprop("test", 1, 0),
        )

    def test_add_material(self):
        coil = MagneticMaterial("coil", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)

        self.assertEqual(
            "mi_addmaterial('coil', 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)",
            FemmWriter().add_material(coil),
        )

        # Electrostatics
        writer = FemmWriter()
        writer.field = femm_electrostatic
        mat = ElectrostaticMaterial("Teflon", 2.1, 2.1, 0)
        self.assertEqual('ei_addmaterial("Teflon", 2.1, 2.1, 0)', writer.add_material(mat))

        # Heat Flow
        writer.field = femm_heat_flow
        mat = HeatFlowMaterial("barack", 1, 2, 3, 4)
        self.assertEqual('hi_addmaterial("barack", 1, 2, 3, 4)', writer.add_material(mat))

        # Current Flow
        writer.field = femm_current_flow
        mat = CurrentFlowMaterial("ribizli", 1, 2, 3, 4, 5, 6)
        self.assertEqual(
            'ci_addmaterial("ribizli", 1, 2, 3, 4, 5, 6)',
            writer.add_material(mat),
        )

    def test_add_boundary(self):
        # mi_addboundprop('abc', 0, 0, 0, 0, 0, 0, 1 / (r * 0.0254 * pi * 4.e-7), 0, 2);

        # dirichlet boundary condition
        dirichlet_boundary = MagneticDirichlet("dirichlet", 1, 2, 3, 4)
        self.assertEqual(
            "mi_addboundprop('dirichlet', 1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0)",
            FemmWriter().add_boundary(dirichlet_boundary),
        )

        # mixed boundary condition
        mixed_boundary = MagneticMixed("mixed_test", 1, 2)
        self.assertEqual(
            "mi_addboundprop('mixed_test', 0, 0, 0, 0, 0, 0, 1, 2, 2, 0, 0)",
            FemmWriter().add_boundary(mixed_boundary),
        )

        # Heatflow tests
        writer = FemmWriter()
        writer.field = femm_heat_flow

        ht_bc = HeatFlowFixedTemperature("Alma", 110)
        self.assertEqual(
            'hi_addboundprop("Alma", 0, 110, 0, 0, 0, 0)',
            writer.add_boundary(ht_bc),
        )

        ht_bc = HeatFlowHeatFlux("Alma", 27)
        self.assertEqual(
            'hi_addboundprop("Alma", 1, 0, 27, 0, 0, 0)',
            writer.add_boundary(ht_bc),
        )

        ht_bc = HeatFlowConvection("Alma", 27, 33)
        self.assertEqual(
            'hi_addboundprop("Alma", 2, 0, 0, 33, 27, 0)',
            writer.add_boundary(ht_bc),
        )

        ht_bc = HeatFlowRadiation("Alma", 66, 27)
        self.assertEqual(
            'hi_addboundprop("Alma", 3, 0, 0, 27, 0, 66)',
            writer.add_boundary(ht_bc),
        )

        ht_bc = HeatFlowPeriodic("Alma")
        self.assertEqual(
            'hi_addboundprop("Alma", 4, 0, 0, 0, 0, 0)',
            writer.add_boundary(ht_bc),
        )

        ht_bc = HeatFlowAntiPeriodic("Alma")
        self.assertEqual(
            'hi_addboundprop("Alma", 5, 0, 0, 0, 0, 0)',
            writer.add_boundary(ht_bc),
        )

        # Electrostatic tests
        writer = FemmWriter()
        writer.field = femm_electrostatic

        el_bc = ElectrostaticFixedVoltage("eper", 10)
        self.assertEqual(
            'ei_addboundprop("eper", 10, 0, 0, 0, 0)',
            writer.add_boundary(el_bc),
        )

        el_bc = ElectrostaticMixed("eper", 1, 9)
        self.assertEqual('ei_addboundprop("eper", 0, 0, 1, 9, 1)', writer.add_boundary(el_bc))

        el_bc = ElectrostaticSurfaceCharge("eper", 156)
        self.assertEqual(
            'ei_addboundprop("eper", 0, 156, 0, 0, 2)',
            writer.add_boundary(el_bc),
        )

        el_bc = ElectrostaticPeriodic("eper")
        self.assertEqual('ei_addboundprop("eper", 0, 0, 0, 0, 3)', writer.add_boundary(el_bc))

        el_bc = ElectrostaticAntiPeriodic("eper")
        self.assertEqual('ei_addboundprop("eper", 0, 0, 0, 0, 4)', writer.add_boundary(el_bc))

        # Current Flow
        writer.field = femm_current_flow
        el_bc = CurrentFlowFixedVoltage("alma", 10)
        self.assertEqual(
            'ci_addboundprop("alma", 10, 0, 0, 0, 0)',
            writer.add_boundary(el_bc),
        )

        el_bc = CurrentFlowMixed("alma", 40, 50)
        self.assertEqual(
            'ci_addboundprop("alma", 0, 0, 40, 50, 2)',
            writer.add_boundary(el_bc),
        )

        el_bc = CurrentFlowSurfaceCurrent("alma", 33)
        self.assertEqual(
            'ci_addboundprop("alma", 0, 33, 0, 0, 2)',
            writer.add_boundary(el_bc),
        )

        el_bc = CurrentFlowPeriodic("alma")
        self.assertEqual('ci_addboundprop("alma", 0, 0, 0, 0, 3)', writer.add_boundary(el_bc))

        el_bc = CurrentFlowAntiPeriodic("alma")
        self.assertEqual('ci_addboundprop("alma", 0, 0, 0, 0, 4)', writer.add_boundary(el_bc))

    def test_addpointprop(self):
        writer = FemmWriter()
        writer.field = femm_electrostatic
        self.assertEqual(
            'ei_addpointprop("alma", 10, 0)',
            writer.add_pointprop("alma", Vp=10),
        )
        self.assertEqual(
            'ei_addpointprop("abc", 0, 0.324)',
            writer.add_pointprop("abc", qp=0.324),
        )

        writer.field = femm_heat_flow
        self.assertEqual(
            'hi_addpointprop("barack", 10, 0)',
            writer.add_pointprop("barack", Tp=10),
        )
        self.assertEqual(
            'hi_addpointprop("cba", 0, 0.324)',
            writer.add_pointprop("cba", qp=0.324),
        )

        writer.field = femm_current_flow
        self.assertEqual(
            'ci_addpointprop("dinnye", 10, 0)',
            writer.add_pointprop("dinnye", Vp=10),
        )
        self.assertEqual(
            'ci_addpointprop("qwert", 0, 0.324)',
            writer.add_pointprop("qwert", qp=0.324),
        )

        writer.field = femm_magnetic
        self.assertEqual(
            'mi_addpointprop("ribizli", 10, 0)',
            writer.add_pointprop("ribizli", a=10),
        )
        self.assertEqual(
            'mi_addpointprop("qwert123", 0, 0.324)',
            writer.add_pointprop("qwert123", j=0.324),
        )
        self.assertEqual(
            'mi_addpointprop("qwert123", 0, 0)',
            writer.add_pointprop("qwert123"),
        )

    def test_block_prop(self):
        self.assertEqual(
            "mi_setblockprop('coil', 0, 0.05, 'icoil', 0, 0, 100)",
            FemmWriter().set_blockprop(
                "coil",
                0,
                0.05,
                0,
                circuit_name="icoil",
                turns=100,
                magdirection=0,
            ),
        )

        writer = FemmWriter()
        writer.field = femm_electrostatic
        self.assertEqual(
            'ei_setblockprop("alma", 1, 2, 3)',
            writer.set_blockprop("alma", 1, 2, 3),
        )

        writer.field = femm_current_flow
        self.assertEqual(
            'ci_setblockprop("alma", 4, 5, 6)',
            writer.set_blockprop("alma", 4, 5, 6),
        )

        writer.field = femm_heat_flow
        self.assertEqual(
            'hi_setblockprop("alma", 7, 8, 9)',
            writer.set_blockprop("alma", 7, 8, 9),
        )

    def test_setarcsegment(self):
        self.assertEqual(
            "mi_setarcsegmentprop(5, 'abc', 0, 0)",
            FemmWriter().set_arc_segment_prop(5, "abc", 0, 0),
        )

    def test_setpointprop(self):
        writer = FemmWriter()
        writer.field = femm_electrostatic
        self.assertEqual('ei_setnodeprop("eper", 0, "<None>")', writer.set_pointprop("eper"))
        self.assertEqual(
            'ei_setnodeprop("eper", 0, "abc")',
            writer.set_pointprop("eper", inductor="abc"),
        )
        self.assertEqual(
            'ei_setnodeprop("eper", 34, "<None>")',
            writer.set_pointprop("eper", groupno=34),
        )

        writer.field = femm_heat_flow
        self.assertEqual(
            'hi_setnodeprop("alma", 23, "barack")',
            writer.set_pointprop("alma", 23, "barack"),
        )
        self.assertEqual('hi_setnodeprop("alma", 0, "<None>")', writer.set_pointprop("alma"))

        writer.field = femm_current_flow
        self.assertEqual('ci_setnodeprop("eper", 0, "<None>")', writer.set_pointprop("eper"))
        self.assertEqual(
            'ci_setnodeprop("eper", 0, "abc")',
            writer.set_pointprop("eper", inductor="abc"),
        )
        self.assertEqual(
            'ci_setnodeprop("eper", 34, "<None>")',
            writer.set_pointprop("eper", groupno=34),
        )

        writer.field = femm_magnetic
        self.assertEqual(
            'mi_setnodeprop("alma", 23, "barack")',
            writer.set_pointprop("alma", 23, "barack"),
        )
        self.assertEqual('mi_setnodeprop("alma", 0, "<None>")', writer.set_pointprop("alma"))

    def test_setsegmentprop(self):
        writer = FemmWriter()
        writer.field = femm_heat_flow
        self.assertEqual(
            'hi_setsegmentprop("alma", 1, 0, 1, 0, "<None>")',
            writer.set_segment_prop("alma", 1, 0, 1, 0),
        )

        writer.field = femm_magnetic
        self.assertEqual(
            'mi_setsegmentprop("eper", 1, 0, 1, 0, "abc")',
            writer.set_segment_prop("eper", 1, 0, 1, 0, "abc"),
        )

        writer.field = femm_electrostatic
        self.assertEqual(
            'ei_setsegmentprop("barak", 0, 0.2, 1, 0, "<None>")',
            writer.set_segment_prop("barak", 0, 0.2, 1, 0),
        )

        writer.field = femm_current_flow
        self.assertEqual(
            'ci_setsegmentprop("ribizli", 0, 0.13, 1, 20, "cba")',
            writer.set_segment_prop("ribizli", 0, 0.13, 1, 20, "cba"),
        )

    def test_run_analysis(self):
        writer = FemmWriter()
        writer.field = femm_electrostatic
        self.assertEqual("ei_analyze(1)", writer.analyze(1))

        writer.field = femm_heat_flow
        self.assertEqual("hi_analyze(0)", writer.analyze(0))

        writer.field = femm_current_flow
        self.assertEqual("ci_analyze(1)", writer.analyze())

        writer.field = femm_magnetic
        self.assertEqual("mi_analyze(2)", writer.analyze(2))

    def test_save_as_command(self):
        writer = FemmWriter()
        writer.field = femm_electrostatic
        self.assertIn("ei_saveas(", writer.save_as("test"))

        writer.field = femm_current_flow
        self.assertIn("ci_saveas(", writer.save_as("test"))

        writer.field = femm_magnetic
        self.assertIn("mi_saveas(", writer.save_as("test"))

        writer.field = femm_heat_flow
        self.assertIn("hi_saveas(", writer.save_as("test"))

    def test_get_circuit_name(self):
        self.assertEqual(
            "current, volt, flux = mo_getcircuitproperties('icoil')",
            FemmWriter().get_circuit_properties("icoil"),
        )

    def test_load_solution(self):
        writer = FemmWriter()
        writer.field = femm_electrostatic
        self.assertEqual("ei_loadsolution()", writer.load_solution())

        writer.field = femm_current_flow
        self.assertEqual("ci_loadsolution()", writer.load_solution())

        writer.field = femm_magnetic
        self.assertEqual("mi_loadsolution()", writer.load_solution())

        writer.field = femm_heat_flow
        self.assertEqual("hi_loadsolution()", writer.load_solution())

    def test_line_integral(self):
        self.assertEqual("mo_lineintegral(0)", FemmWriter().line_integral(0))

    def test_block_integral(self):
        self.assertEqual("mo_blockintegral(30)", FemmWriter().block_integral(30))

    def test_get_point_values(self):
        self.assertEqual("mo_getpointvalues(0.01, 0)", FemmWriter().get_point_values(0.01, 0))

    def test_create_geometry(self):
        """create basic objects: nodes, lines and a circle arc to test the basic functionality of the command."""

        geo = Geometry()

        # test nodes
        a = obj.Node(0.0, 0.0, id=1)
        b = obj.Node(0.0, 1.0, id=2)
        c = obj.Node(1.0, 0.0, id=3)

        geo.nodes = [a, b, c]

        geo.lines = [
            obj.Line(start_pt=a, end_pt=b, id=4),
            obj.Line(start_pt=a, end_pt=c, id=5),
        ]
        geo.circle_arcs = [obj.CircleArc(start_pt=c, center_pt=a, end_pt=b)]

        cmds = FemmWriter().create_geometry(geo)

        self.assertIn("mi_addnode(0.0, 0.0)", cmds)
        self.assertIn("mi_addsegment(0.0, 0.0, 0.0, 1.0)", cmds)
        self.assertIn("mi_addarc(1.0, 0.0, 0.0, 1.0, 90.0, 1)", cmds)

        print(cmds)


# class TestFemmExecutor(TestCase):
#     # def test_executor_proper_file(self):
#     #     testfile = str(Path(__file__).parent / "test.lua")
#     #     warnings.simplefilter("ignore", ResourceWarning)
#     #     exec = FemmExecutor()
#     #     with open(testfile, "w") as f:
#     #         f.write("quit()")
#     #     self.assertEqual(True, exec.run_femm(testfile))

#     def test_executor_invalid_file(self):
#         testfile = str(Path(__file__).parent / "test_invalid.lua")
#         warnings.simplefilter("ignore", ResourceWarning)
#         exec = FemmExecutor()

#         with open(testfile, "w") as f:
#             f.write("not_existing_command()")

#         self.assertEqual(None, exec.run_femm(testfile, timeout=0))
