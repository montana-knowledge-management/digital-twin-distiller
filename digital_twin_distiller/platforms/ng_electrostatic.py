import subprocess
from math import pi
from threading import Timer

import matplotlib.pyplot as plt
import networkx as nx

from digital_twin_distiller import CircleArc, Line, Material, Node
from digital_twin_distiller.boundaries import BoundaryCondition
from digital_twin_distiller.metadata import NgElectrostaticMetadata
from digital_twin_distiller.platforms.platform import Platform
from digital_twin_distiller.utils import get_right_left, pairwise


class NgElectrostatics(Platform):
    def __init__(self, m: NgElectrostaticMetadata):
        super().__init__(m)

        self.G = nx.Graph()
        self.H = nx.DiGraph()

        self.edge_attribures = {"type": None, "rightdomain": 0, "leftdomain": 0, "bc": -1, "mat": -1}

        # materials
        self.mat = {}

    def __copy__(self):
        pass

    def comment(self, str_, nb_newline=1):
        self.file_script_handle.write(f"# {str_}\n")
        self.newline(nb_newline)

    def export_preamble(self):
        ...

    def export_metadata(self):
        ...

    def export_material_definition(self, mat: Material):
        self.mat[mat.name] = mat

    def export_block_label(self, x, y, mat: Material):
        ...

    def export_boundary_definition(self, b: BoundaryCondition):
        ...

    def export_geometry_element(self, e, boundary=None):
        if isinstance(e, Node):
            self.G.add_node(e)

        if isinstance(e, Line):
            attributes = self.edge_attribures.copy()
            attributes["type"] = "line"
            self.G.add_edge(e.start_pt, e.end_pt, **attributes)

        if isinstance(e, CircleArc):
            attributes = self.edge_attribures.copy()
            attributes["type"] = "arc"
            attributes["center_pt"] = tuple(e.center_pt)
            self.G.add_edge(e.start_pt, e.end_pt, **attributes)

    def export_solving_steps(self):
        ...

    def export_results(self, action, entity, variable):
        ...

    def export_closing_steps(self):
        ...

    def execute(self, cleanup=False, timeout=10):
        self.compose_geometry()

        self.open()
        self.ng_export()
        self.close()

        proc = subprocess.Popen(
            ["netgen", self.file_script_handle.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        timer = Timer(timeout, proc.kill)
        try:
            timer.start()
            # stdout, stderr = proc.communicate()
            proc.communicate()
        finally:
            timer.cancel()

    ######################################################################

    def ng_export(self):
        self.ng_export_preamble()
        self.ng_export_metadata()
        self.ng_export_material_definitions()
        self.ng_export_boundary_definitions()
        self.ng_export_geometry()
        self.ng_export_block_labels()
        self.ng_export_solving_steps()
        self.ng_export_postprocessing()
        self.ng_export_closing_steps()

    def ng_export_preamble(self):
        self.comment("PREAMBLE", 1)
        self.write("from ngsolve import *")
        self.write("from netgen.geom2d import SplineGeometry", 2)
        self.write("geo = SplineGeometry()")
        self.newline(2)

    def ng_export_metadata(self):
        self.comment("METADATA", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_material_definitions(self):
        self.comment("MATERIAL DEFINITIONS", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_block_labels(self):
        self.comment("BLOCK LABELS", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_boundary_definitions(self):
        self.comment("BOUNDARY DEFINITIONS", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_geometry(self):
        self.comment("GEOMETRY", 1)

        nodecounter = 0
        nodes = {}

        for ni in self.H.nodes:
            self.write(f"geo.AppendPoint({ni.x}, {ni.y})")
            nodes[ni] = nodecounter
            nodecounter += 1

        for ei, attr in self.H.edges.items():
            ni, nj = ei
            self.write(
                f'geo.Append(["line", {nodes[ni]}, {nodes[nj]}],'
                f' leftdomain={attr["leftdomain"]},'
                f' rightdomain={attr["rightdomain"]},'
                f' bc="gnd")'
            )

        self.write(f"mesh = Mesh(geo.GenerateMesh(maxh=0.1))")

        self.newline(2)

    def ng_export_solving_steps(self):
        self.comment("SOLVER", 1)
        self.write(f"fes = H1(mesh, order=3, dirichlet=[1,2,3,4])")

        self.write("u = fes.TrialFunction()")
        self.write("v = fes.TestFunction()")
        self.write("f = LinearForm(fes)")
        self.write("f += 32 * (y*(1-y)+x*(1-x)) * v * dx")

        self.write("a = BilinearForm(fes, symmetric=True)")
        self.write("a += grad(u)*grad(v)*dx")

        self.write("a.Assemble()")
        self.write("f.Assemble()")

        self.write("gfu = GridFunction(fes)")
        self.write('gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec')

        self.write("Draw (gfu)")

        self.newline(2)

    def ng_export_postprocessing(self):
        self.comment("POSTPROCESSING AND EXPORTING", 1)
        self.comment("empty", 1)

        self.newline(2)

    def ng_export_closing_steps(self):
        self.comment("CLOSING STEPS", 1)
        self.comment("empty", 1)

        self.newline(2)

    ###################################################################

    def render_geo(self):
        plt.figure()
        for nodes, attrs in self.H.edges.items():
            n1, n2 = nodes
            cp = n1.mean(n2)
            ul = 0.5 * n1.unit_to(n2).rotate(pi / 2)
            ur = 0.5 * n1.unit_to(n2).rotate(-pi / 2)

            plt.text(*(cp + ul), attrs["leftdomain"], horizontalalignment="right", verticalalignment="bottom")
            plt.text(*(cp + ur), attrs["rightdomain"], horizontalalignment="left", verticalalignment="top")

            plt.plot([n1.x, n2.x], [n1.y, n2.y], "k-")
            plt.scatter(*n1, c="red", zorder=12)
            plt.scatter(*n2, c="red", zorder=12)

        # render material labels
        for name, mat_i in self.mat.items():
            for xi, yi in mat_i.assigned:
                plt.scatter(xi, yi, c="green", s=80, zorder=100)
                plt.text(xi + 0.1, yi + 0.1, name, verticalalignment="bottom", horizontalalignment="left")

        plt.show()

    def compose_geometry(self):
        material_counter = 1
        mat_i = self.mat["pvc"]
        mat_label = Node(*mat_i.assigned[0])
        start_node = min(self.G.nodes, key=lambda ni: ni.distance_to(mat_label))
        cycles = nx.cycle_basis(self.G, start_node)
        for cycle in cycles:
            edges = pairwise(cycle, includelast=True)
            n0, n1 = next(edges)
            side = get_right_left(n0, n1, mat_label)
            attributes = {**self.G[n0][n1], side: material_counter}
            self.H.add_edge(n0, n1, **attributes)

            for n_start, n_end in edges:
                # __import__('pudb').set_trace()
                attributes = {**self.G[n_start][n_end], side: material_counter}
                self.H.add_edge(n_start, n_end, **attributes)
