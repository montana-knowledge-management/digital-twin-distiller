"""
The goal of this module is to write out the given geometry into a FEMM's lua
script.  The code generated one snapshot from the created model, which can be
run during the Optimization.

The original FEMM code has separate scripting commands for the geometry
generation in different subfields

"""
import subprocess
from collections import namedtuple
from math import asin, degrees
from pathlib import Path
from string import Template
from sys import platform
from threading import Timer

# keywords
femm_current_flow = "current_flow"
femm_electrostatic = "electrostatic"
femm_magnetic = "magnetic"
femm_heat_flow = "heat_flow"

femm_fields = [
    femm_electrostatic,
    femm_magnetic,
    femm_current_flow,
    femm_heat_flow,
]

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
        # Hysteresis lag angle in degrees, used for nonlinear BH curves.
        "Lam_d",
        "Phi_hmax",
        "lam_fill",
        # Fraction of the volume occupied per lamination that is actually
        # filled with iron (Note that this parameter
        # defaults to 1 the femme preprocessor dialog box because, by default,
        # iron completely fills the volume)
        "LamType",
        "Phi_hx",
        # Hysteresis lag in degrees in the x-direction for linear problems.
        "Phi_hy",
        # Hysteresis lag in degrees in the y-direction for linear problems.
        "NStrands",
        # Number of strands in the wire build. Should be 1 for Magnet or Square
        # wire.
        "WireD",
    ],
)  # Diameter of each wire constituent strand in millimeters.

HeatFlowMaterial = namedtuple(
    "heatflow",
    [
        "material_name",
        "kx",  # Thermal conductivity in the x- or r-direction.
        "ky",  # Thermal conductivity in the y- or z-direction.
        "qv",  # Volume heat generation density in units of W/m3.
        "kt",  # Volumetric heat capacity in units of MJ/(m3 * K).
    ],
)

ElectrostaticMaterial = namedtuple(
    "electrostatic",
    [
        "material_name",
        "ex",  # Relative permittivity in the x- or r-direction.
        "ey",  # Relative permittivity in the y- or z-direction.
        "qv",  # Volume charge density in units of C / m3
    ],
)

CurrentFlowMaterial = namedtuple(
    "current_flow",
    [
        "material_name",
        # Electrical conductivity in the x- or r-direction in units of S/m.
        "ox",
        # Electrical conductivity in the y- or z-direction in units of S/m.
        "oy",
        # Relative permittivity in the x- or r-direction.
        "ex",
        # Relative permittivity in the y- or z-direction.
        "ey",
        # Dielectric loss tangent in the x- or r-direction.
        "ltx",
        # Dielectric loss tangent in the y- or z-direction.
        "lty",
    ],
)

# Magnetic Boundary Conditions
MagneticDirichlet = namedtuple("magnetic_dirichlet", ["name", "a_0", "a_1", "a_2", "phi"])
MagneticMixed = namedtuple("magnetic_mixed", ["name", "c0", "c1"])
MagneticAnti = namedtuple("magnetic_anti", ["name"])
MagneticPeriodic = namedtuple("magnetic_periodic", ["name"])
MagneticAntiPeriodicAirgap = namedtuple("magnetic_antiperiodic_airgap", ["name", "angle"])
MagneticPeriodicAirgap = namedtuple("magnetic_periodic_airgap", ["name", "angle"])

# HeatFlow Boundary Conditions
HeatFlowFixedTemperature = namedtuple("heatflow_fixed_temperature", ["name", "Tset"])
HeatFlowHeatFlux = namedtuple("heatflow_heat_flux", ["name", "qs"])
HeatFlowConvection = namedtuple("heatflow_convection", ["name", "h", "Tinf"])
HeatFlowRadiation = namedtuple("heatflow_radiation", ["name", "beta", "Tinf"])
HeatFlowPeriodic = namedtuple("heatflow_periodic", ["name"])
HeatFlowAntiPeriodic = namedtuple("heatflow_anti_periodic", ["name"])

# Electrostatic Boundary Conditions
ElectrostaticFixedVoltage = namedtuple("electrostatic_fixed_voltage", ["name", "Vs"])
ElectrostaticMixed = namedtuple("electrostatic_mixed", ["name", "c0", "c1"])
ElectrostaticSurfaceCharge = namedtuple("electrostatic_surface_charge", ["name", "qs"])
ElectrostaticPeriodic = namedtuple("electrostatic_periodic", ["name"])
ElectrostaticAntiPeriodic = namedtuple("electrostatic_antiperiodic", ["name"])

# Current Flow Boundary Conditions
CurrentFlowFixedVoltage = namedtuple("currentflow_fixed_voltage", ["name", "Vs"])
CurrentFlowMixed = namedtuple("currentflow_mixed", ["name", "c0", "c1"])
CurrentFlowSurfaceCurrent = namedtuple("currentflow_surface_current", ["name", "qs"])
CurrentFlowPeriodic = namedtuple("currentflow_periodidic", ["name"])
CurrentFlowAntiPeriodic = namedtuple("currentflow_antiperiodic", ["name"])


class FemmWriter:
    """Writes out a model snapshot"""

    push = True

    def __init__(self):

        self.field = femm_magnetic
        self.lua_model = []
        self.out_file = "femm_data.csv"

    def validate_field(self, shouldbe=None):
        if self.field not in femm_fields:
            raise ValueError(f"The physical field ({self.field}) is not defined!")

        if shouldbe and shouldbe != self.field:
            raise ValueError(f"({self.field}) != {shouldbe}")

        return True

    def validate_units(self, unit):
        if unit not in {
            "inches",
            "millimeters",
            "centimeters",
            "mils",
            "meters",
            "micrometers",
        }:
            raise ValueError(f"There is no {unit} unit.")
        return True

    def write(self, file_name):
        """Generate a runnable lua-script for a FEMM calculation.

        :param file_name: the code (re)writes the snapshot from the created
                          geometry to the given code
        """

        with open(file_name, "w") as writer:
            for line in self.lua_model:
                writer.write(line + "\n")

    def create_geometry(self, geometry):
        """Creates a FEMM geometry with lua file from the model geometry.

        Building patterns can be:
            - nodes,
            - line segments
            - circle_arcs.

        The field type should be defined separately.
        """

        lua_geometry = []

        # 1 - generate the nodes
        for node in geometry.nodes:
            lua_geometry.append(self.add_node(node.x, node.y))

        for line in geometry.lines:
            lua_geometry.append(
                self.add_segment(
                    line.start_pt.x,
                    line.start_pt.y,
                    line.end_pt.x,
                    line.end_pt.y,
                )
            )

        for arc in geometry.circle_arcs:
            # calculate the angle for the femm circle arc generation
            radius = arc.start_pt.distance_to(arc.center_pt)
            clamp = arc.start_pt.distance_to(arc.end_pt) / 2.0

            deg = 2 * round(degrees(asin(clamp / radius)), 2)

            lua_geometry.append(
                self.add_arc(
                    arc.start_pt.x,
                    arc.start_pt.y,
                    arc.end_pt.x,
                    arc.end_pt.y,
                    angle=deg,
                    maxseg=1,
                )
            )

        return lua_geometry

    def init_problem(self, out_file="femm_data.csv"):
        """
        This commands initialize a femm console and flush the variables
        :param out_file: defines the default output file
        """
        out_file = str(Path(out_file).resolve().as_posix())
        cmd_list = []
        # does nothing if the console is already displayed
        # cmd_list.append("showconsole()")
        # clears both the input and output windows for a fresh start.
        # cmd_list.append("clearconsole()")
        cmd_list.append(f'remove("{out_file}")')  # get rid of the old data file, if it exists
        if self.field == femm_magnetic:
            cmd_list.append("newdocument(0)")  # the 0 specifies a magnetics problem
        if self.field == femm_electrostatic:
            cmd_list.append("newdocument(1)")  # the 1 specifies electrostatics problem
        if self.field == femm_heat_flow:
            cmd_list.append("newdocument(2)")  # the 2 specifies heat flow problem
        if self.field == femm_current_flow:
            cmd_list.append("newdocument(3)")  # the 3 specifies current flow problem

        # cmd_list.append("mi_hidegrid()")
        cmd = Template('file_out = openfile("$outfile", "w")')
        cmd = cmd.substitute(outfile=out_file)
        cmd_list.append(cmd)

        if FemmWriter.push:
            self.lua_model.extend(cmd_list)

        return cmd_list

    def close(self):

        cmd_list = []
        cmd_list.append("closefile(file_out)")
        if self.field == femm_magnetic:
            cmd_list.append("mo_close()")
            cmd_list.append("mi_close()")

        if self.field == femm_heat_flow:
            cmd_list.append("ho_close()")
            cmd_list.append("hi_close()")

        if self.field == femm_electrostatic:
            cmd_list.append("eo_close()")
            cmd_list.append("ei_close()")

        if self.field == femm_current_flow:
            cmd_list.append("co_close()")
            cmd_list.append("ci_close()")

        cmd_list.append("quit()")

        if FemmWriter.push:
            self.lua_model.extend(cmd_list)

        return cmd_list

    def analyze(self, flag=1):
        """
        Runs a FEMM analysis to solve a problem. By default the analysis runs
        in non-visible mode.

        The flag parameter controls whether the fkern window is visible or
        minimized. For a visible window, either specify no value for flag or
        specify 0. For a minimized window, flag should be set to 1.
        """
        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_analyze($flag)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_analyze($flag)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_analyze($flag)")

        if self.field == femm_current_flow:
            cmd = Template("ci_analyze($flag)")

        if FemmWriter.push:
            self.lua_model.append(cmd.substitute(flag=flag))

        return cmd.substitute(flag=flag)

    # object add remove commnads from FEMM MANUAL page 84.
    def add_node(self, x, y):
        """adds a node to the given point (x,y)"""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_addnode($x_coord, $y_coord)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_addnode($x_coord, $y_coord)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_addnode($x_coord, $y_coord)")

        if self.field == femm_current_flow:
            cmd = Template("ci_addnode($x_coord, $y_coord)")

        cmd = cmd.substitute(x_coord=x, y_coord=y)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def add_segment(self, x1, y1, x2, y2, push=True):
        """
        Add a new line segment from node closest to (x1,y1) to node closest to
        (x2,y2)
        """

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_addsegment($x1_coord, $y1_coord, $x2_coord, $y2_coord)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_addsegment($x1_coord, $y1_coord, $x2_coord, $y2_coord)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_addsegment($x1_coord, $y1_coord, $x2_coord, $y2_coord)")

        if self.field == femm_current_flow:
            cmd = Template("ci_addsegment($x1_coord, $y1_coord, $x2_coord, $y2_coord)")

        cmd = cmd.substitute(x1_coord=x1, y1_coord=y1, x2_coord=x2, y2_coord=y2)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def add_blocklabel(self, x, y):
        """Add a new block label at (x,y)"""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_addblocklabel($x_coord, $y_coord)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_addblocklabel($x_coord, $y_coord)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_addblocklabel($x_coord, $y_coord)")

        if self.field == femm_current_flow:
            cmd = Template("ci_addblocklabel($x_coord, $y_coord)")

        cmd = cmd.substitute(x_coord=x, y_coord=y)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def add_arc(self, x1, y1, x2, y2, angle, maxseg):
        """
        Add a new arc segment from the nearest nodeto (x1,y1) to the nearest
        node to (x2,y2) with angle ‘angle’ divided into ‘maxseg’ segments.
        with angle ‘angle’ divided into ‘maxseg’ segments.
        """

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_addarc($x_1, $y_1, $x_2, $y_2, $angle, $maxseg)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_addarc($x_1, $y_1, $x_2, $y_2, $angle, $maxseg)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_addarc($x_1, $y_1, $x_2, $y_2, $angle, $maxseg)")

        if self.field == femm_current_flow:
            cmd = Template("ci_addarc($x_1, $y_1, $x_2, $y_2, $angle, $maxseg)")

        cmd = cmd.substitute(x_1=x1, y_1=y1, x_2=x2, y_2=y2, angle=angle, maxseg=maxseg)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def delete_selected(self):
        """Delete all selected objects"""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = "mi_deleteselected"

        if self.field == femm_electrostatic:
            cmd = "ei_deleteselected"

        if self.field == femm_heat_flow:
            cmd = "hi_deleteselected"

        if self.field == femm_current_flow:
            cmd = "ci_deleteselected"

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def add_boundary(self, boundary):
        """
        :param boundary: checks the type of the boundary parameter, then
        """
        cmd = None
        self.validate_field()

        if self.field == femm_magnetic and isinstance(boundary, MagneticDirichlet):
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

        if self.field == femm_magnetic and isinstance(boundary, MagneticMixed):
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

        if self.field == femm_magnetic and isinstance(boundary, MagneticAnti):
            cmd = Template(
                "mi_addboundprop($propname, $A0, $A1, $A2, $Phi, " "$Mu, $Sig, $c0, $c1, $BdryFormat, $ia, $oa)"
            )
            cmd = cmd.substitute(
                propname="'" + boundary.name + "'",
                A0=0,
                A1=0,
                A2=0,
                Phi=0,
                Mu=0,
                Sig=0,
                c0=0,
                c1=0,
                BdryFormat=5,
                ia=0,
                oa=0,
            )

        if self.field == femm_magnetic and isinstance(boundary, MagneticPeriodic):
            cmd = Template(
                "mi_addboundprop($propname, $A0, $A1, $A2, $Phi, $Mu, $Sig, $c0, $c1, $BdryFormat, $ia, $oa)"
            )
            cmd = cmd.substitute(
                propname="'" + boundary.name + "'",
                A0=0,
                A1=0,
                A2=0,
                Phi=0,
                Mu=0,
                Sig=0,
                c0=0,
                c1=0,
                BdryFormat=4,
                ia=0,
                oa=0,
            )
        if self.field == femm_magnetic and isinstance(boundary, MagneticAntiPeriodicAirgap):
            cmd = Template(
                "mi_addboundprop($propname, $A0, $A1, $A2, $Phi, $Mu, $Sig, $c0, $c1, $BdryFormat, $ia, $oa)"
            )
            cmd = cmd.substitute(
                propname="'" + boundary.name + "'",
                A0=0,
                A1=0,
                A2=0,
                Phi=0,
                Mu=0,
                Sig=0,
                c0=0,
                c1=0,
                BdryFormat=7,
                ia=0,
                oa=boundary.angle,
            )

        if self.field == femm_magnetic and isinstance(boundary, MagneticPeriodicAirgap):
            cmd = Template(
                "mi_addboundprop($propname, $A0, $A1, $A2, $Phi, $Mu, $Sig, $c0, $c1, $BdryFormat, $ia, $oa)"
            )
            cmd = cmd.substitute(
                propname="'" + boundary.name + "'",
                A0=0,
                A1=0,
                A2=0,
                Phi=0,
                Mu=0,
                Sig=0,
                c0=0,
                c1=0,
                BdryFormat=6,
                ia=0,
                oa=boundary.angle,
            )

        # HEATFLOW

        if self.field == femm_heat_flow and isinstance(boundary, HeatFlowFixedTemperature):
            cmd = Template("hi_addboundprop($propname, $BdryFormat, $Tset, $qs, $Tinf, $h, $beta)")
            cmd = cmd.substitute(
                propname=f'"{boundary.name}"',
                BdryFormat=0,
                Tset=boundary.Tset,
                qs=0,
                Tinf=0,
                h=0,
                beta=0,
            )

        if self.field == femm_heat_flow and isinstance(boundary, HeatFlowHeatFlux):
            cmd = Template("hi_addboundprop($propname, $BdryFormat, $Tset, $qs, $Tinf, $h, $beta)")
            cmd = cmd.substitute(
                propname=f'"{boundary.name}"',
                BdryFormat=1,
                Tset=0,
                qs=boundary.qs,
                Tinf=0,
                h=0,
                beta=0,
            )

        if self.field == femm_heat_flow and isinstance(boundary, HeatFlowConvection):
            cmd = Template("hi_addboundprop($propname, $BdryFormat, $Tset, $qs, $Tinf, $h, $beta)")
            cmd = cmd.substitute(
                propname=f'"{boundary.name}"',
                BdryFormat=2,
                Tset=0,
                qs=0,
                Tinf=boundary.Tinf,
                h=boundary.h,
                beta=0,
            )

        if self.field == femm_heat_flow and isinstance(boundary, HeatFlowRadiation):
            cmd = Template("hi_addboundprop($propname, $BdryFormat, $Tset, $qs, $Tinf, $h, $beta)")
            cmd = cmd.substitute(
                propname=f'"{boundary.name}"',
                BdryFormat=3,
                Tset=0,
                qs=0,
                Tinf=boundary.Tinf,
                h=0,
                beta=boundary.beta,
            )

        if self.field == femm_heat_flow and isinstance(boundary, HeatFlowPeriodic):
            cmd = Template("hi_addboundprop($propname, $BdryFormat, $Tset, $qs, $Tinf, $h, $beta)")
            cmd = cmd.substitute(
                propname=f'"{boundary.name}"',
                BdryFormat=4,
                Tset=0,
                qs=0,
                Tinf=0,
                h=0,
                beta=0,
            )

        if self.field == femm_heat_flow and isinstance(boundary, HeatFlowAntiPeriodic):
            cmd = Template("hi_addboundprop($propname, $BdryFormat, $Tset, $qs, $Tinf, $h, $beta)")
            cmd = cmd.substitute(
                propname=f'"{boundary.name}"',
                BdryFormat=5,
                Tset=0,
                qs=0,
                Tinf=0,
                h=0,
                beta=0,
            )

        # Electrostatics
        if self.field == femm_electrostatic and isinstance(boundary, ElectrostaticFixedVoltage):
            cmd = f'ei_addboundprop("{boundary.name}", {boundary.Vs}, 0, 0, 0, 0)'

        if self.field == femm_electrostatic and isinstance(boundary, ElectrostaticMixed):
            cmd = f'ei_addboundprop("{boundary.name}", 0, 0, {boundary.c0}, {boundary.c1}, 1)'

        if self.field == femm_electrostatic and isinstance(boundary, ElectrostaticSurfaceCharge):
            cmd = f'ei_addboundprop("{boundary.name}", 0, {boundary.qs}, 0, 0, 2)'

        if self.field == femm_electrostatic and isinstance(boundary, ElectrostaticPeriodic):
            cmd = f'ei_addboundprop("{boundary.name}", 0, 0, 0, 0, 3)'

        if self.field == femm_electrostatic and isinstance(boundary, ElectrostaticAntiPeriodic):
            cmd = f'ei_addboundprop("{boundary.name}", 0, 0, 0, 0, 4)'

        # Current Flow
        if self.field == femm_current_flow and isinstance(boundary, CurrentFlowFixedVoltage):
            cmd = f'ci_addboundprop("{boundary.name}", {boundary.Vs}, 0, 0, 0, 0)'

        if self.field == femm_current_flow and isinstance(boundary, CurrentFlowMixed):
            cmd = f'ci_addboundprop("{boundary.name}", 0, 0, {boundary.c0}, {boundary.c1}, 2)'

        if self.field == femm_current_flow and isinstance(boundary, CurrentFlowSurfaceCurrent):
            cmd = f'ci_addboundprop("{boundary.name}", 0, {boundary.qs}, 0, 0, 2)'

        if self.field == femm_current_flow and isinstance(boundary, CurrentFlowPeriodic):
            cmd = f'ci_addboundprop("{boundary.name}", 0, 0, 0, 0, 3)'

        if self.field == femm_current_flow and isinstance(boundary, CurrentFlowAntiPeriodic):
            cmd = f'ci_addboundprop("{boundary.name}", 0, 0, 0, 0, 4)'

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def add_material(self, material):
        """
        mi addmaterial("materialname", mu_x, mu_y, H_c, J, Cduct, Lam_d,
                       Phi_hmax, lam_fill, LamType, Phi_hx, Phi_hy,
                       NStrands, WireD)
        """

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
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

        if self.field == femm_electrostatic:
            cmd = Template("ei_addmaterial($materialname, $ex, $ey, $qv)")
            cmd = cmd.substitute(
                materialname=f'"{material.material_name}"',
                ex=material.ex,
                ey=material.ey,
                qv=material.qv,
            )

        if self.field == femm_heat_flow:
            cmd = Template("hi_addmaterial($materialname, $kx, $ky, $qv, $kt)")
            cmd = cmd.substitute(
                materialname=f'"{material.material_name}"',
                kx=material.kx,
                ky=material.ky,
                qv=material.qv,
                kt=material.kt,
            )

        if self.field == femm_current_flow:
            cmd = Template("ci_addmaterial($materialname, $ox, $oy, $ex, $ey, $ltx, $lty)")
            cmd = cmd.substitute(
                materialname=f'"{material.material_name}"',
                ox=material.ox,
                oy=material.oy,
                ex=material.ex,
                ey=material.ey,
                ltx=material.ltx,
                lty=material.lty,
            )

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def add_bhcurve(self, material: MagneticMaterial, file_name, first_column, separator):
        """
        This function adds nonlinear property to an already added linear
        material.

        :param material: The material that is being assigned as nonlinear.
        :param file_name: The location of the file that contains the h-b points
        :param first_column: Assign the first column as either "h" or "b"
        :param separator: The separator that is used in the file
        """
        cmd = ""
        with open(file_name) as f:
            h = 0
            b = 0
            if first_column.lower() == "h":
                h = 0
                b = 1
            else:
                h = 1
                b = 0

            for line in f.readlines():
                try:
                    line = [float(li) for li in line.strip().split(separator)]
                    cmd += f'mi_addbhpoint("{material.material_name}", {line[b]}, {line[h]})\n'
                except ValueError:
                    continue

        if self.push:
            self.lua_model.append(cmd)
        else:
            return cmd

    def delete_selected_nodes(self):
        """
        Delete all selected nodes, the object should be selected the node
        selection command.
        """

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = "mi_deleteselectednodes"

        if self.field == femm_electrostatic:
            cmd = "ei_deleteselectednodes"

        if self.field == femm_heat_flow:
            cmd = "hi_deleteselectednodes"

        if self.field == femm_current_flow:
            cmd = "ci_deleteselectednodes"

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def delete_selected_labels(self):
        """Delete all selected labels"""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = "mi_deleteselectedlabels"

        if self.field == femm_electrostatic:
            cmd = "ei_deleteselectedlabels"

        if self.field == femm_heat_flow:
            cmd = "hi_deleteselectedlabels"

        if self.field == femm_current_flow:
            cmd = "ci_deleteselectedlabels"

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def delete_selected_segments(self):
        """Delete all selected segments."""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = "mi_deleteselectedsegments"

        if self.field == femm_electrostatic:
            cmd = "ei_deleteselectedsegments"

        if self.field == femm_heat_flow:
            cmd = "hi_deleteselectedsegments"

        if self.field == femm_current_flow:
            cmd = "ci_deleteselectedsegments"

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def delete_selected_arc_segments(self):
        """Delete all selected arc segments."""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = "mi_deleteselectedarcsegments"

        if self.field == femm_electrostatic:
            cmd = "ei_deleteselectedarcsegments"

        if self.field == femm_heat_flow:
            cmd = "hi_deleteselectedarcsegments"

        if self.field == femm_current_flow:
            cmd = "ci_deleteselectedarcsegments"

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def add_pointprop(self, propname, **kwargs):
        """
        Adds new point property of name "propname" with various attributes such as:
        Electrostatics:
            - Vp: specified potential Vp (V)
            - qp: point charge density qp (C / m)

        Magnetics:
            - a: specified potential a (Wb / m)
            - j: point current j (A)

        Heat Flow:
            - Tp: specified temperatire Tp at the point (K)
            - qp: point heat generation (W / m)

        Current Flow:
            - Vp: specified potential Vp (V)
            - qp: point current density qp (A / m)
        """
        cmd = None
        # Electrostatics
        if self.field == femm_electrostatic:
            Vp = kwargs.get("Vp", 0)
            qp = kwargs.get("qp", 0)
            cmd = f'ei_addpointprop("{propname}", {Vp}, {qp})'

        # Magnetics
        if self.field == femm_magnetic:
            a = kwargs.get("a", 0)
            j = kwargs.get("j", 0)
            cmd = f'mi_addpointprop("{propname}", {a}, {j})'

        # Heat Flow
        if self.field == femm_heat_flow:
            Tp = kwargs.get("Tp", 0)
            qp = kwargs.get("qp", 0)
            cmd = f'hi_addpointprop("{propname}", {Tp}, {qp})'

        # Current Flow
        if self.field == femm_current_flow:
            Vp = kwargs.get("Vp", 0)
            qp = kwargs.get("qp", 0)
            cmd = f'ci_addpointprop("{propname}", {Vp}, {qp})'

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def add_circprop(self, circuitname, i, circuittype):
        """
        Adds a new circuit property with name "circuitname" with a prescribed
        current, i.  The circuittype parameter is

        Only in the case of magnetic fields.

        :param circuitname: name of the magnetic circuit
        :param i : prescribed current in Amper
        :param circuittype: 0 for a parallel - connected circuit and 1 for a
                            series-connected circuit.
        """
        cmd = f'mi_addcircprop("{circuitname}",{i},{circuittype})'

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    # object selection commnads from FEMM MANUAL page 84.
    def clear_selected(self):
        """Clear all selected nodes, blocks, segments and arc segments."""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = "mi_clearselected()"

        if self.field == femm_electrostatic:
            cmd = "ei_clearselected()"

        if self.field == femm_heat_flow:
            cmd = "hi_clearselected()"

        if self.field == femm_current_flow:
            cmd = "ci_clearselected()"

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def select_segment(self, x, y):
        """Select the line segment closest to (x,y)"""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_selectsegment($xp, $yp)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_selectsegment($xp, $yp)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_selectsegment($xp, $yp)")

        if self.field == femm_current_flow:
            cmd = Template("ci_selectsegment($xp, $yp)")

        cmd = cmd.substitute(xp=x, yp=y)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def select_arc_segment(self, x, y):
        """Select the arc segment closest to (x,y)"""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_selectarcsegment($xp, $yp)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_selectarcsegment($xp, $yp)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_selectarcsegment($xp, $yp)")

        if self.field == femm_current_flow:
            cmd = Template("ci_selectarcsegment($xp, $yp)")

        cmd = cmd.substitute(xp=x, yp=y)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def select_node(self, x, y):
        """Select node closest to (x,y), Returns the coordinates ofthe se-lected node"""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_selectnode($xp, $yp)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_selectnode($xp, $yp)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_selectnode($xp, $yp)")

        if self.field == femm_current_flow:
            cmd = Template("ci_selectnode($xp, $yp)")

        cmd = cmd.substitute(xp=x, yp=y)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def select_label(self, x, y):
        """Select the label closet to (x,y). Returns the coordinates of the selected label."""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_selectlabel($xp, $yp)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_selectlabel($xp, $yp)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_selectlabel($xp, $yp)")

        if self.field == femm_current_flow:
            cmd = Template("ci_selectlabel($xp, $yp)")

        cmd = cmd.substitute(xp=x, yp=y)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def select_group(self, n):
        """
        Select the n th group of nodes, segments, arc segments and block
        labels.  This function will clear all previously selected elements and
        leave the edit mode in 4(group)
        """

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_selectgroup($np)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_selectgroup($np)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_selectgroup($np)")

        if self.field == femm_current_flow:
            cmd = Template("ci_selectgroup($np)")

        cmd = cmd.substitute(np=n)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def select_circle(self, x, y, R, editmode):
        """
        Select circle selects objects within a circle of radius R centered
        at(x, y).If only x, y, and R paramters are given, the current edit mode
        is used.If the editmode parameter is used, 0 denotes nodes, 2 denotes
        block labels, 2 denotes segments, 3 denotes arcs, and 4 specifies that
        all entity types are to be selected.
        """

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_selectcircle($xp, $yp, $Rp, $Editmode)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_selectcircle($xp, $yp, $Rp, $Editmode)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_selectcircle($xp, $yp, $Rp, $Editmode)")

        if self.field == femm_current_flow:
            cmd = Template("ci_selectcircle($xp, $yp, $Rp, $Editmode)")

        cmd = cmd.substitute(xp=x, yp=y, Rp=R, Editmode=editmode)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def select_rectangle(self, x1, y1, x2, y2, editmode):
        """
        This command selects objects within a rectangle definedby points (x1,y1) and (x2,y2).
        If no editmode parameter is supplied, the current edit mode isused. If the editmode parameter is used,
        0 denotes nodes, 2 denotes block labels, 2 denotessegments, 3 denotes arcs, and 4 specifies that all
        entity types are to be selected.
        """

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_selectrectangle($x1p,$y1p,$x2p,$y2p,$Editmode)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_selectrectangle($x1p,$y1p,$x2p,$y2p,$Editmode)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_selectrectangle($x1p,$y1p,$x2p,$y2p,$Editmode)")

        if self.field == femm_current_flow:
            cmd = Template("ci_selectrectangle($x1p,$y1p,$x2p,$y2p,$Editmode)")

        cmd = cmd.substitute(x1p=x1, y1p=y1, x2p=x2, y2p=y2, Editmode=editmode)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def set_pointprop(self, propname, groupno=0, inductor="<None>"):
        """
        :param propname: Set the selected nodes to have the nodal property 'propname'
        :param groupno: Set the selected nodes to have the group number 'groupno'
        :param inductor: Specifies which conductor the node belongs to. Default value is '<None>'
        """
        prefix = None
        cmd = None
        if self.field == femm_magnetic:
            prefix = "mi"
        elif self.field == femm_heat_flow:
            prefix = "hi"
        elif self.field == femm_current_flow:
            prefix = "ci"
        elif self.field == femm_electrostatic:
            prefix = "ei"

        cmd = f'{prefix}_setnodeprop("{propname}", {groupno}, "{inductor}")'

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def set_segment_prop(
        self,
        propname,
        elementsize=1,
        automesh=1,
        hide=0,
        group=0,
        inductor="<None>",
    ):
        """
        :param propname: boundary property
        :param elementsize: Local element size along segment no greater than
                            elementsize.
        :param automesh: mesher defers to the element constraint defined by
                         elementsize, 1 = mesher automatically chooses mesh
                         size along the selected segments
        :param hide: 0 = not hidden in post-processor, 1 == hidden in post
                     processor
        :param group: A member of group number group
        :param inductor: A member of the conductor specified by the string
                         "inconductor". If the segment is not part of a
                         conductor, this parameter can be specified as
                         "<None>".
        """
        prefix = None
        if self.field == femm_magnetic:
            prefix = "mi"
        elif self.field == femm_heat_flow:
            prefix = "hi"
        elif self.field == femm_current_flow:
            prefix = "ci"
        elif self.field == femm_electrostatic:
            prefix = "ei"

        cmd = f'{prefix}_setsegmentprop("{propname}", {elementsize}, {automesh}, {hide}, {group}, "{inductor}")'

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def set_arc_segment_prop(self, maxsegdeg, propname, hide, group):
        """
        :param maxsegdeg: Meshed with elements that span at most maxsegdeg degrees per element
        :param propname: boundary property
        :param hide: 0 = not hidden in post-processor, 1 == hidden in post processor
        :param group: a member of group number group
        """
        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mi_setarcsegmentprop($maxsegdeg, $propname, $hide, $group)")
            cmd = cmd.substitute(
                maxsegdeg=maxsegdeg,
                propname="'" + propname + "'",
                hide=hide,
                group=group,
            )

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def set_blockprop(self, blockname, automesh=1, meshsize=1, group=0, **kwargs):
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
        circuit_name = kwargs.get("circuit_name", "<None>")
        magdirection = kwargs.get("magdirection", 0)
        turns = kwargs.get("turns", 0)

        self.validate_field()

        if self.field == femm_magnetic:
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

        if self.field == femm_heat_flow:
            cmd = Template("hi_setblockprop($blockname, $automesh, $meshsize, $group)")
            cmd = cmd.substitute(
                blockname=f'"{blockname}"',
                automesh=automesh,
                meshsize=meshsize,
                group=group,
            )

        if self.field == femm_electrostatic:
            cmd = Template("ei_setblockprop($blockname, $automesh, $meshsize, $group)")
            cmd = cmd.substitute(
                blockname=f'"{blockname}"',
                automesh=automesh,
                meshsize=meshsize,
                group=group,
            )

        if self.field == femm_current_flow:
            cmd = Template("ci_setblockprop($blockname, $automesh, $meshsize, $group)")
            cmd = cmd.substitute(
                blockname=f'"{blockname}"',
                automesh=automesh,
                meshsize=meshsize,
                group=group,
            )

        if FemmWriter.push:
            self.lua_model.append(cmd)

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

        self.validate_field(femm_magnetic)
        self.validate_units(unit)

        cmd = Template("mi_probdef($frequency,$units,$type,$precision, $depth, $minangle, $acsolver)")
        cmd = cmd.substitute(
            frequency=freq,
            units=r"'" + unit + r"'",
            type=r"'" + type + r"'",
            precision=precision,
            depth=depth,
            minangle=minangle,
            acsolver=acsolver,
        )

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def heat_problem(
        self,
        units,
        type,
        precision=1e-8,
        depth=1,
        minangle=30,
        prevsoln=None,
        timestep=1e-3,
    ):
        """
        :param units: "inches", "millimeters", "centimeters", "mils", "meters", "micrometers"
        :param type: "planar", "axi",
        :param precision: Precision required by the solver. Default value is 1E-8
        :param depth: Depth of the problem into the page for 2D problems
        :param minangle: Minimum angle constraint sen to the mesh generator
        :param prevsoln: Indicates the solution from the previous time step assuming transient time problems
        """
        cmd = None
        self.validate_field(femm_heat_flow)

        self.validate_units(units)

        if type not in {"planar", "axi"}:
            raise ValueError(f"Choose either 'planar' or 'axi', not {type}. ")

        if not prevsoln:
            prevsoln = ""
            timestep = 0

        cmd = f'hi_probdef("{units}", "{type}", {precision}, {depth}, {minangle}, "{prevsoln}", {timestep})'

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def electrostatic_problem(self, units, type, precision=1e-8, depth=1, minangle=30):
        """
        :param units: "inches", "millimeters", "centimeters", "mils", "meters", "micrometers"
        :param type: "planar", "axi",
        :param precision: Precision required by the solver. Default value is 1E-8
        :param depth: Depth of the problem into the page for 2D problems
        :param minangle: Minimum angle constraint sen to the mesh generator
        """

        cmd = None
        self.validate_field(femm_electrostatic)

        self.validate_units(units)

        if type not in {"planar", "axi"}:
            raise ValueError(f"Choose either 'planar' or 'axi', not {type}. ")

        cmd = f'ei_probdef("{units}", "{type}", {precision}, {depth}, {minangle})'

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def currentflow_problem(self, units, type, frequency=0, precision=1e-8, depth=1, minangle=30):
        # TODO: add docstring
        """
        -
        """

        cmd = None
        self.validate_field(femm_current_flow)

        self.validate_units(units)

        if type not in {"planar", "axi"}:
            raise ValueError(f"Choose either 'planar' or 'axi', not {type}. ")

        cmd = f'ci_probdef("{units}", "{type}", {frequency}, {precision}, {depth}, {minangle})'

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def save_as(self, file_name):
        """
        To solve the problem with FEMM, you have to save it with the save_as
        command.

        mi_saveas("filename") saves the file with name "filename". Note if you
        use a path you must use two backslashes e.g. "c:\\temp\\myfemmfile.fem
        """

        cmd = None
        self.validate_field()
        file_name = str(Path(file_name).resolve().as_posix())

        if self.field == femm_magnetic:
            cmd = Template("mi_saveas($filename)")

        if self.field == femm_heat_flow:
            cmd = Template("hi_saveas($filename)")

        if self.field == femm_electrostatic:
            cmd = Template("ei_saveas($filename)")

        if self.field == femm_current_flow:
            cmd = Template("ci_saveas($filename)")

        cmd = cmd.substitute(filename='"' + file_name + '"')

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def load_solution(self):
        """Loads  and displays the solution."""

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = "mi_loadsolution()"

        if self.field == femm_heat_flow:
            cmd = "hi_loadsolution()"

        if self.field == femm_electrostatic:
            cmd = "ei_loadsolution()"

        if self.field == femm_current_flow:
            cmd = "ci_loadsolution()"

        if FemmWriter.push:
            self.lua_model.append(cmd)

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

        self.validate_field()

        if self.field == femm_magnetic:
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

             This function returns one (possibly complex) value,e.g.:volume =
             moblockintegral(10)
        """

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mo_blockintegral($type)")
            cmd = cmd.substitute(type=type)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def get_point_values(self, x, y):
        """
        Get the values associated with the point at x,y Return in order

        Symbol      Definition
        ------      ----------
        A           vector potential A or flux φ
        B1          flux density Bx if planar, Brif axisymmetric
        B2          flux density By if planar, Bzif axisymmetric
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

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("mo_getpointvalues($x, $y)")
            cmd = cmd.substitute(x=x, y=y)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def get_circuit_properties(self, circuit_name, result="current, volt, flux"):
        """Used primarily to obtain impedance information associated with circuit properties.
        Properties are returned for the circuit property named "circuit".
        Three values are returned by the function.

        In order, these results are current, volt and flux of the circuit.
        """

        cmd = None
        self.validate_field()

        if self.field == femm_magnetic:
            cmd = Template("$result = mo_getcircuitproperties($circuit)")

        cmd = cmd.substitute(circuit="'" + circuit_name + "'", result=result)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd

    def write_out_result(self, key, value):
        # writes out a key_value pair
        cmd = Template("write(file_out, '$key', ', ', $value, \"\\{}\")".format("n"))
        cmd = cmd.substitute(key=key, value=value)

        if FemmWriter.push:
            self.lua_model.append(cmd)

        return cmd


class FemmExecutor:
    """
    The goal of this class is to provide a simple and easily configurable FEMM executor.
    This executor uses the Filelink option of FEMM, becuase the slowliness of this file based communication is not critical
    in the case of larger computations, which can be parallelized by ArtapOptimization, or other optimizatino frameworks.
    """

    # Default value of the femm path under linux and under windows.

    home = str(Path.home())
    femm_path_linux = home + "/.wine/drive_c/femm42/bin/femm.exe"
    femm_path_windows = r"C:\femm42\bin\femm.exe"
    executable = femm_path_linux

    def run_femm(self, script_file, timeout=10):
        """This function runs the femm simulation via filelink"""

        cmd_list = []
        script_file = Path(script_file).resolve()
        assert script_file.exists(), f"{script_file} does not exists."

        if platform == "linux":
            FemmExecutor.executable = FemmExecutor.femm_path_linux
            cmd_list.append("wine")
            cmd_list.append(FemmExecutor.executable)
            # script_file = os.popen(f'winepath -w "{script_file}"').read().strip()
            cmd_list.append(f"-lua-script={script_file}")

        elif platform == "win32":
            FemmExecutor.executable = FemmExecutor.femm_path_windows
            cmd_list.append(FemmExecutor.executable)
            cmd_list.append(f"-lua-script={script_file}")
            cmd_list.append("-windowhide")

        proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        timer = Timer(timeout, proc.kill)
        try:
            timer.start()
            stdout, stderr = proc.communicate()
        finally:
            timer.cancel()

        # self.script_file = os.path.basename(script_file)
        # # under linux we are using wine to run FEMM
        # if platform == "linux":
        #     self.femm_command = "wine " + self.femm_path_linux
        #
        #     lua_path = os.path.abspath(script_file)
        #
        #     arg = None
        #     if os.path.isfile(lua_path) and platform == "linux":
        #         arg = '"' + os.popen('winepath -w "' + lua_path + '"').read().strip() + '"'
        #
        #     cmd_string = self.femm_command + f" -lua-script={arg}"
        #
        #     # out = subprocess.run(cmd_string, shell=True, capture_output=True)
        #     proc = Popen(shlex.split(cmd_string), stdout=PIPE, stderr=PIPE)
        #     timer = Timer(timeout, proc.kill)
        #     try:
        #         timer.start()
        #         stdout, stderr = proc.communicate()
        #     finally:
        #         timer.cancel()
        #
        #     if proc.returncode != 0:
        #         err = "Unknown error"
        #         return None
        #     #     if proc.stderr is not None:
        #     #         err = f"Cannot run FEMM.\n\n {proc.stderr}"
        #     #         print(err)
        #
        #     # self.problem.logger.error(err)
        #     # raise RuntimeError(err)
        #
        #     else:
        #         return True
