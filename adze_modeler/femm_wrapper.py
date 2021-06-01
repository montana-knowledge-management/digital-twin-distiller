"""
The goal of this module is to write out the given geometry into a FEMM's lua script.
The code generated one snapshot from the created model, which can be run during the Optimization.

The original FEMM code has separate scripting commands for the geometry generation in different subfields

"""
import os
import subprocess
from collections import namedtuple
from string import Template
from sys import platform

# keywords
kw_current_flow = "current_flow"
kw_electrostatic = "electrostatic"
kw_magnetic = "magnetic"
kw_heat_flow = "heat_flow"

fields = [kw_electrostatic, kw_magnetic, kw_current_flow, kw_heat_flow]

# material types for the different FEMM suppoerted magnetic fields

# Lam_type
# 0 – Not laminated or laminated in plane
# 1 – laminated x or r
# 2 – laminated y or z
# 3 – Magnet wire
# 4 – Plain stranded wire
# 5 – Litz wire
# 6 – Square wire

MagneticMaterial = namedtuple(
    "magnetic",
    [
        "material_name",
        "mu_x",  # Relative permeability in the x- or r-direction.
        "mu_y",  # Relative permeability in the y- or z-direction.
        "H_c",  # Permanent magnet coercivity in Amps/Meter.
        "J",  # Real Applied source current density in Amps/mm2 .
        "Cduct",  # Electrical conductivity of the material in MS/m.
        "Lam_d",  # Hysteresis lag angle in degrees, used for nonlinear BH curves.
        "Phi_hmax",
        "lam_fill",
        # Fraction of the volume occupied per lamination that is actually filled with iron (Note that this parameter defaults to 1 the femme preprocessor dialog box because, by default, iron completely fills the volume)
        "LamType",
        "Phi_hx",
        # Hysteresis lag in degrees in the x-direction for linear problems.
        "Phi_hy",
        # Hysteresis lag in degrees in the y-direction for linear problems.
        "NStrands",
        # Number of strands in the wire build. Should be 1 for Magnet or Square wire.
        "WireD",
    ],
)  # Diameter of each wire constituent strand in millimeters.

# TODO: other types should be defined
MagneticDirichlet = namedtuple("magnetic_dirichlet", ["name", "a_0", "a_1", "a_2", "phi"])
MagneticMixed = namedtuple("magnetic_mixed", ["name", "c0", "c1"])


class FemmWriter:
    """Writes out a model snapshot"""

    field = kw_magnetic
    lua_model = []  # list of the lua commands

    def write(self, file_name):
        """Generate a runnable lua-script for a FEMM calculation.

        :param file_name: the code (re)writes the snapshot from the created geometry to the given code
        """

        with open(file_name, "w") as writer:
            for line in self.lua_model:
                writer.write(line + "\n")

    def init_problem(self, out_file="femm_data.csv"):
        """
        This commands initialize a femm console and flush the variables
        :param out_file: defines the default output file
        """

        cmd_list = []
        cmd_list.append("showconsole()")  # does nothing if the console is already displayed
        cmd_list.append("clearconsole()")  # clears both the input and output windows for a fresh start.
        cmd_list.append(f'remove("{out_file}")')  # get rid of the old data file, if it exists
        if self.field == kw_magnetic:
            cmd_list.append("newdocument(0)")  # the 0 specifies a magnetics problem
        elif self.field == kw_electrostatic:
            cmd_list.append("newdocument(1)")  # the 1 specifies electrostatics problem
        elif self.field == kw_heat_flow:
            cmd_list.append("newdocument(2)")  # the 2 specifies heat flow problem
        elif self.field == kw_current_flow:
            cmd_list.append("newdocument(3)")  # the 3 specifies current flow problem


        # cmd_list.append("mi_hidegrid()")
        cmd = Template('file_out = openfile("$outfile", "w")')
        cmd = cmd.substitute(outfile=out_file)
        cmd_list.append(cmd)
        return cmd_list

    def close(self, out_file="femm_data.csv"):

        cmd_list = []
        cmd_list.append("closefile(file_out)")
        cmd_list.append("mo_close()")
        cmd_list.append("mi_close()")
        cmd_list.append("quit()")

        return cmd_list

    def analyze(self, flag=1):
        """

        Runs a FEMM analysis to solve a problem. By default the analysis runs in non-visible mode.

        The flag parameter controls whether the fkern window is visible or minimized. For a visible window,
        either specify no value for flag or specify 0. For a minimized window, flag should be set to 1.
        """
        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_analyze($flag)")

        return cmd.substitute(flag=flag)

    # object add remove commnads from FEMM MANUAL page 84.
    def add_node(self, x, y):
        """adds a node to the given point (x,y)"""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_addnode($x_coord, $y_coord)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_addnode($x_coord, $y_coord)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_addnode($x_coord, $y_coord)")

        if self.field == kw_current_flow:
            cmd = Template("ci_addnode($x_coord, $y_coord)")

        return cmd.substitute(x_coord=x, y_coord=y)

    def add_segment(self, x1, y1, x2, y2):
        """Add a new line segment from node closest to (x1,y1) to node closest to (x2,y2)"""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_addsegment($x1_coord, $y1_coord, $x2_coord, $y2_coord)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_addsegment($x1_coord, $y1_coord, $x2_coord, $y2_coord)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_addsegment($x1_coord, $y1_coord, $x2_coord, $y2_coord)")

        if self.field == kw_current_flow:
            cmd = Template("ci_addsegment($x1_coord, $y1_coord, $x2_coord, $y2_coord)")

        return cmd.substitute(x1_coord=x1, y1_coord=y1, x2_coord=x2, y2_coord=y2)

    def add_blocklabel(self, x, y):
        """Add a new block label at (x,y)"""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_addblocklabel($x_coord, $y_coord)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_addblocklabel($x_coord, $y_coord)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_addblocklabel($x_coord, $y_coord)")

        if self.field == kw_current_flow:
            cmd = Template("ci_addblocklabel($x_coord, $y_coord)")

        return cmd.substitute(x_coord=x, y_coord=y)

    def add_arc(self, x1, y1, x2, y2, angle, maxseg):
        """
        Add a new arc segment from the nearest nodeto (x1,y1) to the nearest node to (x2,y2)
        with angle ‘angle’ divided into ‘maxseg’ segments.
        with angle ‘angle’ divided into ‘maxseg’ segments.
        """

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_addarc($x_1, $y_1, $x_2, $y_2, $angle, $maxseg)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_addarc($x_1, $y_1, $x_2, $y_2, $angle, $maxseg)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_addarc($x_1, $y_1, $x_2, $y_2, $angle, $maxseg)")

        if self.field == kw_current_flow:
            cmd = Template("ci_addarc($x_1, $y_1, $x_2, $y_2, $angle, $maxseg)")

        return cmd.substitute(x_1=x1, y_1=y1, x_2=x2, y_2=y2, angle=angle, maxseg=maxseg)

    def delete_selected(self):
        """Delete all selected objects"""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = "mi_deleteselected"

        if self.field == kw_electrostatic:
            cmd = "ei_deleteselected"

        if self.field == kw_heat_flow:
            cmd = "hi_deleteselected"

        if self.field == kw_current_flow:
            cmd = "ci_deleteselected"

        return cmd

    def add_boundary(self, boundary):
        """
        :param boundary: checks the type of the boundary parameter, then
        """
        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic and isinstance(boundary, MagneticDirichlet):
            cmd = Template(
                "mi_addboundprop($propname, $A0, $A1, $A2, $Phi, $Mu, $Sig, " "$c0, $c1, $BdryFormat, $ia, $oa)"
            )
            cmd = cmd.substitute(
                propname="'" + boundary.name + "'",
                A0=boundary.a_0,
                A1=boundary.a_1,
                A2=boundary.a_2,
                Phi=boundary.phi,
                Mu=0,
                Sig=0,
                c0=0,
                c1=0,
                BdryFormat=0,
                ia=0,
                oa=0,
            )

        if self.field == kw_magnetic and isinstance(boundary, MagneticMixed):
            cmd = Template(
                "mi_addboundprop($propname, $A0, $A1, $A2, $Phi, $Mu, $Sig, " "$c0, $c1, $BdryFormat, $ia, $oa)"
            )
            cmd = cmd.substitute(
                propname="'" + boundary.name + "'",
                A0=0,
                A1=0,
                A2=0,
                Phi=0,
                Mu=0,
                Sig=0,
                c0=boundary.c0,
                c1=boundary.c1,
                BdryFormat=2,
                ia=0,
                oa=0,
            )

        return cmd

    def add_material(self, material):
        """
        mi addmaterial("materialname", mu_x, mu_y, H_c, J, Cduct, Lam_d, Phi_hmax,
                        lam_fill, LamType, Phi_hx, Phi_hy, NStrands, WireD)
        """

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic and isinstance(material, MagneticMaterial):
            cmd = Template(
                "mi_addmaterial($materialname, $mux, $muy, $Hc, $J, $Cduct, $Lamd, $Phi_hmax, $lamfill, "
                "$LamType, $Phi_hx, $Phi_hy, $NStrands, $WireD)"
            )

            # cmd.substitute(x_1=x1, y_1=y1, x_2=x2, y_2=y2, angle=angle, maxseg=maxseg)
            # hy is missing from the FEMM command
            cmd = cmd.substitute(
                materialname="'" + material.material_name + "'",
                mux=material.mu_x,
                muy=material.mu_y,
                Hc=material.H_c,
                J=material.J,
                Cduct=material.Cduct,
                Lamd=material.Lam_d,
                Phi_hmax=material.Phi_hmax,
                lamfill=material.lam_fill,
                LamType=material.LamType,
                Phi_hx=material.Phi_hx,
                Phi_hy=material.Phi_hy,
                NStrands=material.NStrands,
                WireD=material.WireD,
            )
        #
        # if self.field == kw_electrostatic:
        #     pass
        #
        # if self.field == kw_heat_flow:
        #     pass
        #
        # if self.field == kw_current_flow:
        #     pass

        return cmd

    def delete_selected_nodes(self):
        """Delete all selected nodes, the object should be selected the node selection command."""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = "mi_deleteselectednodes"

        if self.field == kw_electrostatic:
            cmd = "ei_deleteselectednodes"

        if self.field == kw_heat_flow:
            cmd = "hi_deleteselectednodes"

        if self.field == kw_current_flow:
            cmd = "ci_deleteselectednodes"

        return cmd

    def delete_selected_labels(self):
        """Delete all selected labels"""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = "mi_deleteselectedlabels"

        if self.field == kw_electrostatic:
            cmd = "ei_deleteselectedlabels"

        if self.field == kw_heat_flow:
            cmd = "hi_deleteselectedlabels"

        if self.field == kw_current_flow:
            cmd = "ci_deleteselectedlabels"

        return cmd

    def delete_selected_segments(self):
        """Delete all selected segments."""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = "mi_deleteselectedsegments"

        if self.field == kw_electrostatic:
            cmd = "ei_deleteselectedsegments"

        if self.field == kw_heat_flow:
            cmd = "hi_deleteselectedsegments"

        if self.field == kw_current_flow:
            cmd = "ci_deleteselectedsegments"

        return cmd

    def delete_selected_arc_segments(self):
        """Delete all selected arc segments."""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = "mi_deleteselectedarcsegments"

        if self.field == kw_electrostatic:
            cmd = "ei_deleteselectedarcsegments"

        if self.field == kw_heat_flow:
            cmd = "hi_deleteselectedarcsegments"

        if self.field == kw_current_flow:
            cmd = "ci_deleteselectedarcsegments"

        return cmd

    def add_circprop(self, circuitname, i, circuittype):
        """
        Adds a new circuit property with name "circuitname" with a prescribed current, i.
        The circuittype parameter is

        Only in the case of magnetic fields.

        :param circuitname: name of the magnetic circuit
        :param i : prescribed current in Amper
        :param circuittype:  0for a parallel - connected circuit and 1 for a series-connected circuit.
        """

        return f'mi_addcircprop("{circuitname}",{i},{circuittype})'

    # object selection commnads from FEMM MANUAL page 84.
    def clear_selected(self):
        """Clear all selected nodes, blocks, segments and arc segments."""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = "mi_clearselected()"

        if self.field == kw_electrostatic:
            cmd = "ei_clearselected()"

        if self.field == kw_heat_flow:
            cmd = "hi_clearselected()"

        if self.field == kw_current_flow:
            cmd = "ci_clearselected()"

        return cmd

    def select_segment(self, x, y):
        """Select the line segment closest to (x,y)"""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_selectsegment($xp, $yp)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_selectsegment($xp, $yp)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_selectsegment($xp, $yp)")

        if self.field == kw_current_flow:
            cmd = Template("ci_selectsegment($xp, $yp)")

        return cmd.substitute(xp=x, yp=y)

    def select_arc_segment(self, x, y):
        """Select the arc segment closest to (x,y)"""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_selectarcsegment($xp, $yp)")

        # if self.field == kw_electrostatic:
        #     cmd = Template("ei_selectnode($xp, $yp)")
        #
        # if self.field == kw_heat_flow:
        #     cmd = Template("hi_selectnode($xp, $yp)")
        #
        # if self.field == kw_current_flow:
        #     cmd = Template("ci_selectnode($xp, $yp)")

        return cmd.substitute(xp=x, yp=y)

    def select_node(self, x, y):
        """Select node closest to (x,y), Returns the coordinates ofthe se-lected node"""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_selectnode($xp, $yp)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_selectnode($xp, $yp)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_selectnode($xp, $yp)")

        if self.field == kw_current_flow:
            cmd = Template("ci_selectnode($xp, $yp)")

        return cmd.substitute(xp=x, yp=y)

    def select_label(self, x, y):
        """Select the label closet to (x,y). Returns the coordinates of the selected label."""

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_selectlabel($xp, $yp)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_selectlabel($xp, $yp)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_selectlabel($xp, $yp)")

        if self.field == kw_current_flow:
            cmd = Template("ci_selectlabel($xp, $yp)")

        return cmd.substitute(xp=x, yp=y)

    def select_group(self, n):
        """
        Select the n th group of nodes, segments, arc segments and block labels.
        This function will clear all previously selected elements and leave the edit mode in 4(group)
        """

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_selectgroup($np)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_selectgroup($np)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_selectgroup($np)")

        if self.field == kw_current_flow:
            cmd = Template("ci_selectgroup($np)")

        return cmd.substitute(np=n)

    def select_circle(self, x, y, R, editmode):
        """
        Select circle selects objects within a circle of radius R centered at(x, y).If only x, y, and R paramters
        are given, the current edit mode is used.If the editmode parameter is used, 0 denotes nodes, 2 denotes block
        labels, 2 denotes segments, 3 denotes arcs, and 4 specifies that all entity types are to be selected.
        """

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_selectcircle($xp, $yp, $Rp, $Editmode)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_selectcircle($xp, $yp, $Rp, $Editmode)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_selectcircle($xp, $yp, $Rp, $Editmode)")

        if self.field == kw_current_flow:
            cmd = Template("ci_selectcircle($xp, $yp, $Rp, $Editmode)")

        return cmd.substitute(xp=x, yp=y, Rp=R, Editmode=editmode)

    def select_rectangle(self, x1, y1, x2, y2, editmode):
        """
        This command selects objects within a rectangle definedby points (x1,y1) and (x2,y2).
        If no editmode parameter is supplied, the current edit mode isused. If the editmode parameter is used,
        0 denotes nodes, 2 denotes block labels, 2 denotessegments, 3 denotes arcs, and 4 specifies that all
        entity types are to be selected.
        """

        cmd = None
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_selectrectangle($x1p,$y1p,$x2p,$y2p,$Editmode)")

        if self.field == kw_electrostatic:
            cmd = Template("ei_selectrectangle($x1p,$y1p,$x2p,$y2p,$Editmode)")

        if self.field == kw_heat_flow:
            cmd = Template("hi_selectrectangle($x1p,$y1p,$x2p,$y2p,$Editmode)")

        if self.field == kw_current_flow:
            cmd = Template("ci_selectrectangle($x1p,$y1p,$x2p,$y2p,$Editmode)")

        return cmd.substitute(x1p=x1, y1p=y1, x2p=x2, y2p=y2, Editmode=editmode)

    def set_arc_segment_prop(self, maxsegdeg, propname, hide, group):
        """
        :param maxsegdeg: Meshed with elements that span at most maxsegdeg degrees per element
        :param propname: boundary property
        :param hide: 0 = not hidden in post-processor, 1 == hidden in post processor
        :param group: a member of group number group
        """
        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_setarcsegmentprop($maxsegdeg, $propname, $hide, $group)")
            cmd = cmd.substitute(maxsegdeg=maxsegdeg, propname="'" + propname + "'", hide=hide, group=group)
        return cmd

    def set_blockprop(self, blockname=None, meshsize=None, circuit_name=None, magdirection=0, group="group", turns=0):
        """
        :param meshsize: default value is None -> invokes automesh
            this command will use automesh option as the default, if the mesh size is not defined

        # these parameters used only in the case of magnetic field

        :param magdirection:

            The magnetization is directed along an angle in measured in degrees denoted by the
            parameter magdirection. Alternatively, magdirection can be a string containing a
            formula that prescribes the magnetization direction as a function of element position.
            In this formula theta and R denotes the angle in degrees of a line connecting the center
            each element with the origin and the length of this line, respectively; x and y denote
            the x- and y-position of the center of the each element. For axisymmetric problems, r
            and z should be used in place of x and y.

        :param group: None, mebmer of the named group

        """
        cmd = None

        if not meshsize:
            automesh = 1
            meshsize = 0

        else:
            automesh = 0

        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template(
                "mi_setblockprop($blockname, $automesh, $meshsize, $incircuit, $magdirection, $group, $turns)"
            )
            cmd = cmd.substitute(
                blockname="'" + blockname + "'",
                automesh=automesh,
                meshsize=meshsize,
                incircuit="'" + circuit_name + "'",
                magdirection=magdirection,
                group=group,
                turns=turns,
            )

        # ei_setblockprop("blockname", automesh, meshsize, group)
        # ci_setblockprop("blockname", automesh, meshsize, group)
        # hi_setblockprop("blockname", automesh, meshsize, group)

        return cmd

    # problem commands for the magnetic problem
    def magnetic_problem(self, freq, unit, type, precision=1e-8, depth=1, minangle=30, acsolver=0):
        """
         Definition of the magnetic problem, like probdef(0,'inches','axi',1e-8,0,30);

         :param freq: Frequency in Hertz (required)
         :param unit: "inches","millimeters","centimeters","mils","meters, and"micrometers" (required)
         :param type: "planar", "axi" (required)
         :param precision: 1e-8 (required)
         :param depth: depth of the analysis (not mandatory)
         :param minangle: sent to the mesh generator to define the minimum angle of the meshing triangles(not mandatory)
         :param acsolver: the selected acsolver for the problem (not mandatory) - 0 successive approximation, 1 Newton solver


        The generated lua command has the following role:

         miprobdef(frequency,units,type,precision,(depth),(minangle),(acsolver) changes the problem definition.
         Set frequency to the desired frequency in Hertz. The units parameter specifies the units used for measuring
         length in the problem domain. Valid"units"en-tries are"inches","millimeters","centimeters","mils","meters,
         and"micrometers".Set the parameter problem type to"planar"for a 2-D planar problem, or to"axi"for
         anaxisymmetric problem. The precision parameter dictates the precision required by the solver.
         For example, entering 1E-8 requires the RMS of the residual to be less than 10−8.A fifth parameter,
         representing the depth of sthe problem in the into-the-page direction for2-D planar problems, can also also be
         specified. A sixth parameter represents the minimumangle constraint sent to the mesh generator, 30 degress is
         the usual choice. The acsolver parameter specifies which solver is to be used for AC problems:
         0 for successive approximation, 1 for Newton. A seventh parameter specifies the solver type tobe used
         for AC problems.
        """

        if self.field != kw_magnetic:
            raise ValueError("Set the magnetic field parameter!")

        cmd = Template("mi_probdef($frequency,$units,$type,$precision, $depth, $minangle, $acsolver)")
        return cmd.substitute(
            frequency=freq,
            units=r"'" + unit + r"'",
            type=r"'" + type + r"'",
            precision=precision,
            depth=depth,
            minangle=minangle,
            acsolver=acsolver,
        )

    def heat_problem(self, units, type, precision=1e-8, depth=1, minangle=30, prevsoln=0, timestep=1e-3):
        if self.field != kw_heat_flow:
            raise ValueError("Set the heat flow problem type!")

        if units not in {"inches", "millimeters", "centimeters", "mils", "meters", "micrometers"}:
            raise ValueError(f"There is no {units} unit.")

        if type not in {"planar", "axi"}:
            raise ValueError(f"Choose either 'planar' or 'axi', not {type}. ")

        return f'hi_probdef("{units}", "{type}", {precision}, {depth}, {minangle}, {prevsoln}, {timestep})'




    def save_as(self, file_name):
        """
        To solve the problem with FEMM, you have to save it with the save_as command.

        mi_saveas("filename") saves the file with name "filename". Note if you use a path you
                              must use two backslashes e.g. "c:\\temp\\myfemmfile.fem
        """

        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mi_saveas($filename)")

        return cmd.substitute(filename='"' + file_name + '"')

    def load_solution(self):
        """Loads  and displays the solution."""

        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = "mi_loadsolution()"

        return cmd

    # post processing commands --- data extraction
    def line_integral(self, type):
        """
        Calculate the line integral of the defined contour.

        Parameter name, values 1, values 2, values 3, values 4
        0 -- Bn ---  total Bn, avg Bn, - , -
        1 -- Ht ---  total Ht, avg Ht, -, -
        2 -- Contour length  --- surface area, - , - , -
        3 -- Stress Tensor Force --- DC r/x force, DC y/z force, 2x r/x force, 2x y/z force
        4 -- Stress Tensor Torque --- total (B.n)^2, avg (B.n)^2
        """

        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mo_lineintegral($type)")
            cmd = cmd.substitute(type=type)

        return cmd

    def block_integral(self, type):
        """
        Calculate the block integral of the selected blocks.

             Type        Definition
            ------      ------------
             0            AJ
             1            A
             2            Magnetic Field Energy
             3            Hysteresis and/or lamination losses
             4            Resistive losses
             5            Block cross-section area
             6            Total losses
             7            Total current
             8            Integral of Bx (Br) over the block
             9            Integral of By (Bax) over the block
             10           Block volume
             11           x (or r) part of steady state Lorentz force
             12           y (or z) part of stead state Lorentz force
             13           x (or r) part of steady state 2x Lorentz force
             14           y (or z) part of stead state 2x Lorentz force
             15           steady state lorentz torque
             16           2×component of Lorentz torque
             17           Magnetic field coenergy
             18           x (or r) part of steady-state weighted stress tensor force
             19           y (or z) part of steady-state weighted stress tensor force
             20           x (or r) part of 2×weighted stress tensor force
             21           y (or z) part of 2×weighted stress tensor force
             22           Steady-state weighted stress tensor torque
             23           2×component of weighted stress tensor torque
             24           R2(i.e.moment of inertia / density)
             25           x (or r) part of 1×weighted stress tensor force
             26           y (or z) part of 1×weighted stress tensor force
             27           1×component of weighted stress tensor torque
             28           x (or r) part of 1×Lorentz force
             29           y (or z) part of 1×Lorentz force
             30           1×component of Lorentz torque

             This function returns one (possibly complex) value,e.g.:volume = moblockintegral(10)
        """

        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mo_blockintegral($type)")
            cmd = cmd.substitute(type=type)

        return cmd

    def get_point_values(self, x, y):
        """
        Get the values associated with the point at x,y Return in order

        Symbol      Definition
        ------      ----------
        A           vector potential A or flux φ
        B1          flux density    Bxif planar,Brif axisymmetric
        B2          flux density Byif planar, Bzif axisymmetric
        Sig         electrical conductivity σ
        E           stored energy density
        H1          field intensity Hxif planar,Hrif axisymmetric
        H2          field intensity Hyif planar,Hzif axisymmetric
        Je          eddy current density
        Js          source current density
        Mu1         relative permeability μxif planar,μrif axisymmetric
        Mu2relative permeability μyif planar,μzif axisymmetric
        Pe          Power density dissipated through ohmic losses
        Ph          Power density dissipated by hysteresis
        """

        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("mo_getpointvalues($x, $y)")
            cmd = cmd.substitute(x=x, y=y)

        return cmd

    def get_circuit_properties(self, circuit_name, result="current, volt, flux"):
        """Used primarily to obtain impedance information associated with circuit properties.
        Properties are returned for the circuit property named "circuit".
        Three values are returned by the function.

        In order, these results are current, volt and flux of the circuit."""

        if self.field not in fields:
            raise ValueError("The physical field is not defined!")

        if self.field == kw_magnetic:
            cmd = Template("$result = mo_getcircuitproperties($circuit)")

        return cmd.substitute(circuit="'" + circuit_name + "'", result=result)

    def write_out_result(self, key, value):
        # writes out a key_value pair
        cmd = Template("write(file_out, '$key', ', ', $value, \"\\{}\") \n".format("n"))
        return cmd.substitute(key=key, value=value)


class FemmExecutor:
    """
    The goal of this class is to provide a simple and easily configurable FEMM executor.
    This executor uses the Filelink option of FEMM, becuase the slowliness of this file based communication is not critical
    in the case of larger computations, which can be parallelized by Artap, or other optimizatino frameworks.
    """

    # Default value of the femm path under linux and under windows.
    femm_path_linux = "$HOME/.wine/drive_c/femm42/bin/femm.exe"
    femm_path_windows = r"C:\FEMM42\bin\femm.exe"

    def run_femm(self, script_file):
        """This function runs the femm simulation via filelink"""

        self.script_file = os.path.basename(script_file)

        # under linux we are using wine to run FEMM
        if platform == "linux":
            self.femm_command = "wine " + self.femm_path_linux

            lua_path = os.path.abspath(self.script_file)

            arg = None
            if os.path.isfile(lua_path) and platform == "linux":
                arg = '"' + os.popen('winepath -w "' + lua_path + '"').read().strip() + '"'

            cmd_string = self.femm_command + f" -lua-script={arg}"

            out = subprocess.run(cmd_string, shell=True, stdout=subprocess.PIPE)

            if out.returncode != 0:
                err = "Unknown error"
                if out.stderr is not None:
                    err = f"Cannot run FEMM.\n\n {out.stderr}"
                    print(err)
                # self.problem.logger.error(err)
                # raise RuntimeError(err)
