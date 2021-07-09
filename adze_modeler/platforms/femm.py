import os
from collections import Iterable
from math import asin
from copy import copy

from adze_modeler.boundaries import BoundaryCondition, DirichletBoundaryCondition
from adze_modeler.material import Material
from adze_modeler.metadata import Metadata
from adze_modeler.objects import Node, Line, CircleArc
from adze_modeler.platforms.platform import Platform
from adze_modeler.femm_wrapper import FemmWriter, FemmExecutor
from adze_modeler.femm_wrapper import  femm_magnetic, femm_electrostatic, femm_heat_flow, femm_current_flow
from adze_modeler.femm_wrapper import MagneticDirichlet
from adze_modeler.femm_wrapper import MagneticMaterial
from glob import glob
from adze_modeler.boundaries import AntiPeriodicBoundaryCondition
from adze_modeler.femm_wrapper import MagneticAnti


class Femm(Platform):

    def __init__(self, m: Metadata):
        super().__init__(m)
        self.writer = FemmWriter()
        self.writer.push = False

        if self.metadata.problem_type == 'magnetic':
            self.writer.field = femm_magnetic

        elif self.metadata.problem_type == 'electrostatic':
            self.writer.field = femm_electrostatic

        elif self.metadata.problem_type == 'heat':
            self.writer.field = femm_heat_flow

        elif self.metadata.problem_type == 'current':
            self.writer.field = femm_current_flow
        else:
            raise ValueError()

    def __copy__(self):
        newplatform = Femm(copy(self.metadata))
        newplatform.writer.field = self.writer.field
        return newplatform

    def comment(self, str_, nb_newline=1):
        self.file_script_handle.write(f"-- {str_}\n")
        self.newline(nb_newline)

    def export_preamble(self):
        pass

    def export_metadata(self):
        cmdlist = self.writer.init_problem(out_file=self.metadata.file_metrics_name)
        for cmd_i in cmdlist:
            self.write(cmd_i)

        type_ = 'axi' if self.metadata.coordinate_type == 'axisymmetric' else 'planar'

        if self.metadata.problem_type == 'magnetic':
            self.write(self.writer.magnetic_problem(freq=self.metadata.frequency,
                                                    unit=self.metadata.unit,
                                                    type=type_,
                                                    precision=self.metadata.precision,
                                                    depth=self.metadata.depth,
                                                    minangle=self.metadata.minangle,
                                                    acsolver=self.metadata.acsolver
                                                    ))

        if not self.metadata.smartmesh:
            self.write("smartmesh(0)")

    def export_material_definition(self, mat: Material):
        if self.metadata.problem_type == 'magnetic':
            femm_material = MagneticMaterial(material_name=mat.name,
                                             mu_x=mat.mu_r,
                                             mu_y=mat.mu_r,
                                             H_c=mat.coercivity,
                                             J=mat.Je / 1.0e6, # A / m2 -> MA / m2
                                             Cduct=mat.conductivity / 1.0e6,
                                             Lam_d=mat.thickness,
                                             lam_fill=mat.fill_factor,
                                             NStrands=0.0,
                                             WireD=0.0,
                                             LamType=0,
                                             Phi_hmax=0,
                                             Phi_hx=0,
                                             Phi_hy=0)

            self.write(self.writer.add_material(femm_material))
            if isinstance(mat.b, Iterable):
                assert len(mat.b) == len(mat.h), "B and H should have the same length"
                for bi, hi in zip(mat.b, mat.h):
                    self.write(f'mi_addbhpoint("{mat.name}", {bi}, {hi})')

    def export_block_label(self, x, y, mat: Material):
        self.write(self.writer.add_blocklabel(x, y))
        self.write(self.writer.select_label(x, y))
        self.write(self.writer.set_blockprop(blockname=mat.name, automesh=int(self.metadata.smartmesh), meshsize=mat.meshsize))
        self.write(self.writer.clear_selected())

    def export_boundary_definition(self, b: BoundaryCondition):
        if isinstance(b, DirichletBoundaryCondition):
            if self.metadata.problem_type == 'magnetic':
                femm_boundary = MagneticDirichlet(name=b.name,
                                                  a_0=b.valuedict['magnetic_potential'],
                                                  a_1=b.valuedict['magnetic_potential'],
                                                  a_2=b.valuedict['magnetic_potential'],
                                                  phi=0
                                                  )

        if isinstance(b, AntiPeriodicBoundaryCondition):
            if self.metadata.problem_type == 'magnetic':
                femm_boundary = MagneticAnti(name=b.name,
                                             )

        self.write(self.writer.add_boundary(femm_boundary))

    def export_geometry_element(self, e, boundary=None):
        automesh = 1
        elementsize = 1
        if self.metadata.elementsize:
            automesh = 0
            elementsize = 1
        else:
            automesh = 1
            elementsize = self.metadata.elementsize

        if isinstance(e, Node):
            self.write(self.writer.add_node(e.x, e.y))

        if isinstance(e, Line):
            self.write(self.writer.add_segment(e.start_pt.x, e.start_pt.y, e.end_pt.x, e.end_pt.y))
            if boundary:
                # we should give an internal point to select the line
                m_x = (e.start_pt.x + e.end_pt.x) * 0.5
                m_y = (e.start_pt.y + e.end_pt.y) * 0.5

                self.write(self.writer.select_segment(m_x, m_y))
                self.write(self.writer.set_segment_prop(boundary, automesh=automesh, elementsize=elementsize))
                self.write(self.writer.clear_selected())

        if isinstance(e, CircleArc):
            # we should find an internal point of the circle arc
            # to achieve this the start node rotated with deg/2

            radius = e.start_pt.distance_to(e.center_pt)
            clamp = e.start_pt.distance_to(e.end_pt) / 2.0
            theta = round(asin(clamp / radius), 2)
            internal_pt = e.start_pt.rotate_about(e.center_pt, theta)

            self.write(self.writer.select_arc_segment(m_x, m_y))
            self.write(self.writer.set_arc_segment_prop(e.max_seg_deg, boundary))
            self.write(self.writer.clear_selected())

    def export_solving_steps(self):
        femm_filename = self.get_script_name()

        if self.metadata.problem_type == 'magnetic':
            femm_filename += '.fem'
        elif self.metadata.problem_type == 'electrostatic':
            femm_filename += '.fee'
        elif self.metadata.problem_type == 'current':
            femm_filename += '.fec'
        elif self.metadata.problem_type == 'heat':
            femm_filename += '.feh'

        self.write(self.writer.save_as(femm_filename))
        self.write(self.writer.analyze())
        self.write(self.writer.load_solution())


    def export_metrics(self, action, entity, variable):
        mappings = {
            "Bx": 'B1',
            "By": 'B2',
            "Br": 'B1',
            "Bz": 'B2',
            "Hx": 'H1',
            "Hy": 'H2',
            "Hr": 'H1',
            "Hz": 'H2',
        }
        if action == 'point_value':
            x = entity[0]
            y = entity[1]
            self.write('A, B1, B2, Sig, E, H1, H2, Je, Js, Mu1, Mu2, Pe, Ph = ', nb_newline=0)
            self.write(self.writer.get_point_values(x, y))
            self.write(f'write(file_out, "{variable}, {x}, {y}, ", {mappings[variable]}, "\\n")')

        if action == "mesh_info":
            fieldmapping = {'electrostatic': 'eo', 'magnetic':'mo', 'heat':'ho', 'current':'co'}
            self.write(f'write(file_out, "nodes, ", {fieldmapping[self.metadata.problem_type]}_numnodes(), "\\n")')
            self.write(f'write(file_out, "elements, ", {fieldmapping[self.metadata.problem_type]}_numelements(), "\\n")')

    def export_closing_steps(self):
        for cmd_i in self.writer.close():
            self.write(cmd_i)
        pass

    def execute(self, cleanup=False):
        executor = FemmExecutor()
        try:
            executor.run_femm(self.metadata.file_script_name)
            if cleanup:
                femm_files = glob(f'{self.get_script_name()}.*')
                for file_i in femm_files:
                    os.remove(file_i)

            return True

        except Exception as e:
            print(e)
            return False
