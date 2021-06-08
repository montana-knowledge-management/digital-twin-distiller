from math import asin

from adze_modeler.femm_wrapper import femm_fields
from adze_modeler.femm_wrapper import femm_magnetic
from adze_modeler.femm_wrapper import FemmWriter
from adze_modeler.objects import CircleArc
from adze_modeler.objects import Line


class BoundarySnaphot:
    """
    The goal of this class is to create and bind the selected boundaries with the exising ones.

    """

    def __init__(self):
        self.field_type = None
        self.boundary = None
        self.elements = []  # This list contains the elements of which defined by the given boundary condition
        self.elementsize = None

        # optional femm commands
        self.max_seg_deg = 1
        self.group = 0
        self.in_circuit = "<None>"

    def create_femm_boundaries(self):
        cmd = []

        # adds the selected Femm command for the defined boundary
        cmd.append(FemmWriter().add_boundary(self.boundary))

        for elem in self.elements:

            if isinstance(elem, Line):
                # we should give an internal point to select the line
                m_x = (elem.start_pt.x + elem.end_pt.x) * 0.5
                m_y = (elem.start_pt.y + elem.end_pt.y) * 0.5

                cmd.append(FemmWriter().select_segment(m_x, m_y))

                if self.elementsize:
                    automesh = 0
                    elementsize = 1
                else:
                    automesh = 1
                    elementsize = self.elementsize

                cmd.append(
                    FemmWriter().set_segment_prop(
                        self.boundary.name,
                        elementsize=elementsize,
                        automesh=automesh,
                        hide=0,
                        group=self.group,
                        inductor=self.in_circuit,
                    )
                )
                cmd.append(FemmWriter().clear_selected())

            if isinstance(elem, CircleArc):
                # we should find an internal point of the circle arc
                # to achieve this the start node rotated with deg/2

                radius = elem.start_pt.distance_to(elem.center_pt)
                clamp = elem.start_pt.distance_to(elem.end_pt) / 2.0

                theta = round(asin(clamp / radius), 2)

                internal_pt = elem.start_pt.rotate_about(elem.center_pt, theta)

                cmd.append(FemmWriter().select_arc_segment(internal_pt.x, internal_pt.y))
                cmd.append(FemmWriter().set_arc_segment_prop(self.max_seg_deg, self.boundary.name, 0, self.group))
                cmd.append(FemmWriter().clear_selected())

        return cmd
