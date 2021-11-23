import subprocess
from threading import Timer

import matplotlib.pyplot as plt
import networkx as nx

from digital_twin_distiller import CircleArc, Line, Material, Node
from digital_twin_distiller.boundaries import BoundaryCondition
from digital_twin_distiller.metadata import NgElectrostaticMetadata
from digital_twin_distiller.platforms.platform import Platform
from digital_twin_distiller.utils import pairwise
from shapely.geometry import Point, LinearRing, Polygon
import itertools as it



class NgElectrostatics(Platform):
    def __init__(self, m: NgElectrostaticMetadata):
        super().__init__(m)

        self.G = nx.Graph()
        self.H = nx.DiGraph()

        self.edge_attribures = {"type": None, "rightdomain": 0, "leftdomain": 0, "bc": -1}

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
        # self.render_geo()
        # exit(0)

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
        self.write(f'fes = H1(mesh, order=3, dirichlet="gnd")')

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
        pos = {ni: (ni.x, ni.y) for ni in self.G.nodes}
        plt.figure(figsize=(10, 10))
        nx.draw(self.H, pos)

        edge_labels = {}
        for u, v, data in self.H.edges(data=True):
            edge_labels[u, v] = f"l - {data['leftdomain']}\nr - {data['rightdomain']}"
        # for li in labels:
        #     plt.scatter(*li, c='k')
        #     plt.text(li.x + 0.05, li.y + 0.05, li.label, horizontalalignment='left')
        nx.draw_networkx_edge_labels(self.H, pos, edge_labels, rotate=False)
        plt.show()

    def compose_geometry(self):
        """
        1. Find all cycles in the undirected graph
        2. Select the cycles that contains only 1 label
        3. Transform the cycles into directed graphs
        4. merge these graphs
        """
        # material_counter = 1
        # mat_i = self.mat["pvc"]
        # mat_label = Node(*mat_i.assigned[0])

        # convert the cycles into a list of undirected graphs
        allcycle = self._find_all_cycles(self.G)

        loops, labels = self._filter_cycles(allcycle)

        surfaces = self._generate_surfaces(loops, labels)

        self.H = nx.DiGraph()
        for si in surfaces:
            for u, v, attr in si.edges(data=True):
                if self.H.has_edge(u, v):
                    attr2 = self.H.get_edge_data(u, v)
                    if attr2['leftdomain'] == 0:
                        self.H[u][v]['leftdomain'] = attr['leftdomain']

                    if attr2['rightdomain'] == 0:
                        self.H[u][v]['rightdomain'] = attr['rightdomain']


                elif self.H.has_edge(v, u):
                    attr2 = self.H.get_edge_data(v, u)
                    if attr2['leftdomain'] == 0:
                        self.H[v][u]['leftdomain'] = attr['r']

                    if attr2['rightdomain'] == 0:
                        self.H[v][u]['rightdomain'] = attr['leftdomain']

                else:
                    self.H.add_edge(u, v, **attr)

    def _find_all_cycles(self, G):
        """
        forked from networkx dfs_edges function. Assumes nodes are integers, or at least
        types which work with min() and > .

        Credit: https://gist.github.com/joe-jordan/6548029
        """

        nodes = [list(i)[0] for i in nx.connected_components(G)]

        # extra variables for cycle detection:
        cycle_stack = []
        output_cycles = set()

        def get_hashable_cycle(cycle):
            """cycle as a tuple in a deterministic order."""
            m = min(cycle)
            mi = cycle.index(m)
            mi_plus_1 = mi + 1 if mi < len(cycle) - 1 else 0
            if cycle[mi - 1] > cycle[mi_plus_1]:
                result = cycle[mi:] + cycle[:mi]
            else:
                result = list(reversed(cycle[:mi_plus_1])) + list(reversed(cycle[mi_plus_1:]))
            return tuple(result)

        for start in nodes:
            if start in cycle_stack:
                continue
            cycle_stack.append(start)

            stack = [(start, iter(G[start]))]
            while stack:
                parent, children = stack[-1]
                try:
                    child = next(children)

                    if child not in cycle_stack:
                        cycle_stack.append(child)
                        stack.append((child, iter(G[child])))
                    else:
                        i = cycle_stack.index(child)
                        if i < len(cycle_stack) - 2:
                            output_cycles.add(get_hashable_cycle(cycle_stack[i:]))

                except StopIteration:
                    stack.pop()
                    cycle_stack.pop()

        allcycle = (tuple(i) for i in output_cycles)

        # convert cycles into graphs
        return tuple(nx.Graph(pairwise(cycle, cycle=True)) for cycle in allcycle)

    def _filter_cycles(self, cycles):
        labelcounter = []
        for Hi in cycles:
            labels_i = []
            lr = Polygon(LinearRing([tuple(ni) for ni in Hi.nodes]))
            for mat_name, mat_i in self.mat.items():
                for lable_position in mat_i.assigned:
                    if lr.contains(Point(lable_position)):
                        labels_i.append(mat_name)

            labelcounter.append(labels_i)

        selector = tuple(map(lambda li: len(li) == 1, labelcounter))
        loops = it.compress(cycles, selector)
        selected_labels = tuple(li[0] for li in it.compress(labelcounter, selector))

        return loops, selected_labels

    def _generate_surfaces(self, loops, labels):
        surfaces = []
        for loop, label in zip(loops, labels):
            si = nx.DiGraph()
            area = 0.0
            edges = nx.find_cycle(loop)
            for n_start, n_end in edges:
                attrs = {'leftdomain': 0, 'rightdomain': 0}
                area += (n_end.x - n_start.x) * (n_end.y + n_start.y)
                si.add_edge(n_start, n_end, **attrs)

            orientation = 'counter-clockwise' if area < 0 else 'clockwise'
            for n_start, n_end in edges:
                if orientation == 'clockwise':
                    si[n_start][n_end]['rightdomain'] = label
                else:
                    si[n_start][n_end]['leftdomain'] = label

            surfaces.append(si)

        return surfaces