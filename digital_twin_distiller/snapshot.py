from collections import defaultdict

from digital_twin_distiller.boundaries import BoundaryCondition
from digital_twin_distiller.geometry import Geometry
from digital_twin_distiller.material import Material
from digital_twin_distiller.platforms.platform import Platform
from digital_twin_distiller.utils import getID


class Snapshot:
    def __init__(self, p: Platform):
        self.id = getID()
        self.platform = p

        self.boundaries = {}
        self.materials = {}
        self.nodes = {}
        self.lines = {}
        self.circle_arcs = {}
        self.metrics = []

    def set_platform(self, p: Platform):
        self.platform = p

    def add_boundary_condition(self, bc: BoundaryCondition):
        if bc.name in self.boundaries.keys():
            raise ValueError("This boundary is already added")
        elif bc.field != self.platform.metadata.problem_type:
            raise TypeError(f"Boundary condition field type != problem field type")
        else:
            self.boundaries[bc.name] = bc

    def assign_boundary_condition(self, x, y, name):

        if name not in self.boundaries.keys():
            raise ValueError(f'There is no boundary condition called "{name}"')

        closest_line = min(self.lines.values(), key=lambda li: li.distance_to_point(x, y))

        self.boundaries[name].assigned.add(closest_line.id)

    def assign_arc_boundary_condition(self, x, y, name):

        if name not in self.boundaries.keys():
            raise ValueError(f'There is no boundary condition called "{name}"')

        closest_arc = min(
            self.circle_arcs.values(),
            key=lambda arc_i: arc_i.distance_to_point(x, y),
        )

        self.boundaries[name].assigned.add(closest_arc.id)

    def add_material(self, mat: Material):
        if mat.name in self.materials.keys():
            raise ValueError("This material is already added")
        else:
            self.materials[mat.name] = mat

    def assign_material(self, x, y, name):
        if name in self.materials.keys():
            self.materials[name].assigned.append((x, y))
        else:
            raise ValueError(f'There is no material called "{name}"')

    def add_geometry(self, geo: Geometry):
        for ni in geo.nodes:
            self.nodes[ni.id] = ni

        for li in geo.lines:
            self.lines[li.id] = li

        for arc_i in geo.circle_arcs:
            self.circle_arcs[arc_i.id] = arc_i

    def add_postprocessing(self, action, entity, variable):
        self.metrics.append((action, entity, variable))

    def export(self, customfilehandle=None, develmode=None):
        self.platform.open(customfilehandle)
        self.platform.export_preamble()

        self.platform.newline(1)
        self.platform.comment("PROBLEM")
        self.platform.export_metadata()

        self.platform.newline(1)
        self.platform.comment("MATERIAL DEFINITIONS")
        for name, mat in self.materials.items():
            self.platform.export_material_definition(mat)

        self.platform.newline(1)
        self.platform.comment("BOUNDARY DEFINITIONS")
        for name, bi in self.boundaries.items():
            self.platform.export_boundary_definition(bi)

        self.platform.newline(1)
        self.platform.comment("GEOMETRY")
        # nodes
        for id, node_i in self.nodes.items():
            self.platform.export_geometry_element(node_i)

        # Export the boundaries first
        for name_i, boundary_i in self.boundaries.items():
            for id_i in boundary_i.assigned:
                if id_i in self.lines.keys():
                    line_i = self.lines.pop(id_i)
                    self.platform.export_geometry_element(line_i, boundary=name_i)
                elif id_i in self.circle_arcs.keys():
                    arc_i = self.circle_arcs.pop(id_i)
                    self.platform.export_geometry_element(arc_i, boundary=name_i)
                else:
                    raise ValueError(
                        "There is no line with the id:",
                        id_i,
                        "boundary:",
                        boundary_i.name,
                    )

        # Export the rest
        for id, line_i in self.lines.items():
            self.platform.export_geometry_element(line_i)

        # Export circle arcs
        for id, arc_i in self.circle_arcs.items():
            self.platform.export_geometry_element(arc_i)

        self.platform.newline(1)
        self.platform.comment("BLOCK LABELS")
        for name_i, material_i in self.materials.items():
            for xi, yi in material_i.assigned:
                self.platform.export_block_label(xi, yi, material_i)

        if not develmode:
            self.platform.newline(1)
            self.platform.comment("SOLVE")
            self.platform.export_solving_steps()

            self.platform.newline(1)
            self.platform.comment("POSTPROCESSING AND EXPORTING")
            for step in self.metrics:
                self.platform.export_metrics(step[0], step[1], step[2])

            self.platform.newline(1)
            self.platform.comment("CLOSING STEPS")
            self.platform.export_closing_steps()

        self.platform.close()

    def __call__(self, *args, **kwargs):
        try:
            self.export()
            # self.execute()
            # return self.retrive_results()
            return 1
        except Exception as e:
            return None

    def execute(self, cleanup=False, timeout=10):
        return self.platform.execute(cleanup=cleanup, timeout=timeout)

    def retrive_results(self):
        results = defaultdict(list)
        with open(self.platform.metadata.file_metrics_name) as f:
            for line in f.readlines():
                line = line.strip().split(",")
                if len(line) == 2:
                    results[line[0]] = float(line[1])
                else:
                    results[line.pop(0)].append(tuple(float(xi) for xi in line))

        return results
