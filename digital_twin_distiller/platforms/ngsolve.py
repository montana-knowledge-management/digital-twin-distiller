import itertools as it
import subprocess
from abc import ABCMeta, abstractmethod
from threading import Timer

import matplotlib.pyplot as plt
import networkx as nx
from shapely.geometry import LinearRing, Point, Polygon

from digital_twin_distiller import CircleArc, Line, Material, Node
from digital_twin_distiller.boundaries import BoundaryCondition
from digital_twin_distiller.metadata import NgSolveMetadata
from digital_twin_distiller.platforms.platform import Platform
from digital_twin_distiller.utils import pairwise


class NgSolve(Platform, metaclass=ABCMeta):
    """
    Realize an abstract class for the ngsolve based platforms. NgSolve lets to create a specific fem solver for your
    custom partial differential equation system.

    Every solution should be defined as a separate platform, like ng_electrostatic.
    """

    def __init__(self, m: NgSolveMetadata):
        super().__init__(m)

        self.G = nx.Graph()
        self.H = nx.DiGraph()

        self.edge_attribures = {"type": None, "rightdomain": 0, "leftdomain": 0, "bc": -1}

        # materials
        self.mat = {}

    def comment(self, str_, nb_newline=1):
        self.file_script_handle.write(f"# {str_}\n")
        self.newline(nb_newline)

    def export_preamble(self):
        """This function is not used here."""
        ...

    def export_metadata(self):
        """This function is not used here."""
        ...

    def export_material_definition(self, mat: Material):
        self.mat[mat.name] = mat

    def export_block_label(self, x, y, mat: Material):
        """This function is not used here."""
        ...

    def export_boundary_definition(self, b: BoundaryCondition):
        """This function is not used here."""
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
        """This function is not used here."""
        ...

    def export_results(self, action, entity, variable, custom_name):
        """This function is not used here."""
        ...

    def export_closing_steps(self):
        """This function is not used here."""
        ...

    def execute(self, cleanup=False, timeout=10, **kwargs):
        self.compose_geometry()
        # self.render_geo()
        # exit(0)

        self.open()
        self.ng_export()
        self.close()

        runner = kwargs.get("runner", "netgen")

        proc = subprocess.Popen([runner, self.file_script_handle.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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

    @abstractmethod
    def ng_export_metadata(self):
        ...

    @abstractmethod
    def ng_export_material_definitions(self):
        ...

    @abstractmethod
    def ng_export_block_labels(self):
        ...

    @abstractmethod
    def ng_export_boundary_definitions(self):
        ...

    def ng_export_geometry(self):
        self.comment("GEOMETRY", 1)

        nodecounter = it.count()
        nodes = {}

        materials = dict(zip(self.mat.keys(), it.count(1)))
        materials[0] = 0

        for ni in self.H.nodes:
            self.write(f"geo.AppendPoint({ni.x}, {ni.y})")
            nodes[ni] = next(nodecounter)

        for ni, nj, attr in self.H.edges(data=True):
            left_material = attr["leftdomain"]
            right_material = attr["rightdomain"]
            self.write(
                f'geo.Append(["line", {nodes[ni]}, {nodes[nj]}],'
                f" leftdomain={materials[left_material]},"
                f" rightdomain={materials[right_material]},"
                f' bc="gnd")'
            )

        self.write(f"mesh = Mesh(geo.GenerateMesh(maxh=0.1))")

        self.newline(2)

    @abstractmethod
    def ng_export_solving_steps(self):
        ...

    @abstractmethod
    def ng_export_postprocessing(self):
        ...

    @abstractmethod
    def ng_export_closing_steps(self):
        ...

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
        # plt.show()

        return plt.gca()

    def compose_geometry(self):
        """
        1. Find all cycles in the undirected graph
        2. Select the cycles that contains only 1 label
        3. Transform the cycles into directed graphs
        4. merge these graphs
        """

        # convert the cycles into a list of undirected graphs
        allcycle = self._find_all_cycles(self.G)

        # filter the cycles
        loops, labels = self._filter_cycles(allcycle)

        # generate surfaces
        surfaces = self._generate_surfaces(loops, labels)

        # merge surfaces into the H directed graph
        self.H = self._merge_surfaces(surfaces)

        return surfaces

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

    def _filter_cycles(self, cycles: list):
        """
        This function filters the cycles based on how many labels one cycle contains.
        """
        labelcounter = []
        for Hi in cycles:
            labels_i = []
            # convert the cycle into a polygon
            lr = Polygon(LinearRing([tuple(ni) for ni in Hi.nodes]))

            # iterate over the materials
            for mat_name, mat_i in self.mat.items():
                # iterate over the materials label points
                for label_position in mat_i.assigned:
                    # if the label is in the polygon append it to the label set
                    if lr.contains(Point(label_position)):
                        labels_i.append(mat_name)

            # append all labels that is in 1 cycle to the labelcounter list
            labelcounter.append(labels_i)

        # labelselector is now a list of lists, now we have to select the elements with length=1
        # selector is now a list of booleans
        selector = tuple(map(lambda li: len(li) == 1, labelcounter))

        # using selector, filter the cycles and labels
        loops = it.compress(cycles, selector)
        selected_labels = tuple(li[0] for li in it.compress(labelcounter, selector))

        return loops, selected_labels

    def _generate_surfaces(self, loops, labels):
        """
        This function turns a list of undirected graphs into a list of
        directed graphs with edge attributes. These edge attributes are the
        leftdomain and rightdomain with the proper material names.
        """
        surfaces = []
        for loop, label in zip(loops, labels):
            # the new surface
            si = nx.DiGraph()

            # only 1 cycle exist
            edges = nx.find_cycle(loop)

            # the signed area will help to decide the orientation of the cycle found above.
            area = 0.0
            for n_start, n_end in edges:
                # copy the default edge attributes
                attrs = self.edge_attribures.copy()

                # add the new edges' contribution to the area
                area += (n_end.x - n_start.x) * (n_end.y + n_start.y)

                # add the edge to the surface
                si.add_edge(n_start, n_end, **attrs)

            # iterate over the edges and set the materials name on the left/right side based on the cycles' orientation
            orientation = "counter-clockwise" if area < 0 else "clockwise"
            for n_start, n_end in edges:
                if orientation == "clockwise":
                    si[n_start][n_end]["rightdomain"] = label
                else:
                    si[n_start][n_end]["leftdomain"] = label

            surfaces.append(si)

        return surfaces

    def _merge_surfaces(self, surfaces: list):
        """
        This function takes a list of directed graph and merges them into 1 directed graph.
        It will detect the edges that are already in the graph and update the edge attributes
        accordingly. This happens regardless of the edges' direction: (u, v) == (v, u)

        """

        T = nx.DiGraph()

        # iterate over the surfaces
        for si in surfaces:

            # iterate over the surface edges
            for u, v, attr in si.edges(data=True):
                if T.has_edge(u, v):
                    attr2 = T.get_edge_data(u, v)
                    if attr2["leftdomain"] == 0:
                        T[u][v]["leftdomain"] = attr["leftdomain"]

                    if attr2["rightdomain"] == 0:
                        T[u][v]["rightdomain"] = attr["rightdomain"]

                elif T.has_edge(v, u):
                    attr2 = T.get_edge_data(v, u)
                    if attr2["leftdomain"] == 0:
                        T[v][u]["leftdomain"] = attr["r"]

                    if attr2["rightdomain"] == 0:
                        T[v][u]["rightdomain"] = attr["leftdomain"]

                else:
                    T.add_edge(u, v, **attr)

        return T
