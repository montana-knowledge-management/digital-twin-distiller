from copy import copy
from numpy import linspace

from adze_modeler.model import BaseModel
from adze_modeler.metadata import FemmMetadata
from adze_modeler.platforms.femm import Femm
from adze_modeler.snapshot import Snapshot
from adze_modeler.material import Material
from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.modelpaths import ModelDir
from adze_modeler.modelpiece import ModelPiece
from adze_modeler.objects import Node, Line, ParametricBezier
from adze_modeler.utils import pairwise

ModelDir.set_base(__file__)


class MagneticClutch(BaseModel):
    """docstring for MagneticClutch"""

    def __init__(self, **kwargs):
        super(MagneticClutch, self).__init__(**kwargs)
        self._init_directories()

    def setup_solver(self):
        femm_metadata = FemmMetadata()
        femm_metadata.problem_type = "magnetic"
        femm_metadata.coordinate_type = "planar"
        femm_metadata.file_script_name = self.file_solver_script
        femm_metadata.file_metrics_name = self.file_solution
        femm_metadata.unit = "millimeters"
        femm_metadata.smartmesh = True
        femm_metadata.depth = 1000

        self.platform = Femm(femm_metadata)
        self.snapshot = Snapshot(self.platform)

    def define_materials(self):
        air = Material('air')

        magnet = Material('N52')
        magnet.mu_r = 1.05
        magnet.coercivity = 956180
        magnet.conductivity = 0.667e6
        magnet.b = [0.0, 0.0732, 0.1464, 0.2195, 0.2927, 0.3659, 0.4391, 1.4636]
        magnet.h = [0.0, 7484, 18446, 37084, 72766, 124571, 179757, 956180]

        magnet_up = copy(magnet)
        magnet_up.name = "N52UP"
        magnet_up.remanence_angle = 90

        magnet_down = copy(magnet)
        magnet_down.name = "N52DOWN"
        magnet_down.remanence_angle = -90

        printed_steel = Material('printed_steel')
        printed_steel.conductivity = 10.44e6
        printed_steel.b = [0.000000, 0.040456, 0.050456, 0.060456, 0.070456, 0.080456, 0.090456, 0.100456, 0.110456,
                           0.120456, 0.130456, 0.140456, 0.150456, 0.160456, 0.170456, 0.180456, 0.190456, 0.200456,
                           0.210456, 0.220456, 0.230456, 0.240456, 0.250456, 0.260456, 0.270456, 0.280456, 0.290456,
                           0.300456, 0.310456, 0.320456, 0.330456, 0.340456, 0.350456, 0.360456, 0.370456, 0.380456,
                           0.390456, 0.400456, 0.410456, 0.420456, 0.430456, 0.440456, 0.450456, 0.460456, 0.470456,
                           0.480456, 0.490456, 0.500456, 0.510456, 0.520456, 0.530456, 0.540456, 0.550456, 0.560456,
                           0.570456, 0.580456, 0.590456, 0.600456, 0.610456, 0.620456, 0.630456, 0.640456, 0.650456,
                           0.660456, 0.670456, 0.680456, 0.690456, 0.700456, 0.710456, 0.720456, 0.730456, 0.740456,
                           0.750456, 0.760456, 0.770456, 0.780456, 0.790456, 0.800456, 0.810456, 0.820456, 0.830456,
                           0.840456, 0.850456, 0.860456, 0.870456, 0.880456, 0.890456, 0.900456, 0.910456, 0.920456,
                           0.930456, 0.940456, 0.950456, 0.960456, 0.970456, 0.980456, 0.990456, 1.000456, 1.010456,
                           1.020456, 1.030456, 1.040456, 1.050456, 1.060456, 1.070456, 1.080456, 1.090456, 1.100456,
                           1.110456, 1.120456, 1.130456, 1.140456, 1.150456, 1.160456, 1.170456, 1.180456, 1.190456,
                           1.200456, 1.210456, 1.288064, 1.298064, 1.308064, 1.318064, 1.328064, 1.338064, 1.348064,
                           1.358064, 1.368064, 1.378064, 1.388064, 1.398064, 1.408064, 1.418064, 1.428064, 1.438064,
                           1.448064, 1.458064, 1.468064, 1.478064, 1.488064, 1.498064, 1.508064, 1.518064, 1.528064,
                           1.538064, 1.548064, 1.558064, 1.568064]
        printed_steel.h = [0.000000, 1.248605, 2.538744, 3.803118, 5.053133, 6.305102, 7.574851, 8.878180, 10.230892,
                           11.648786, 13.147664, 14.738680, 16.413652, 18.140858, 19.882087, 21.599123, 23.255076,
                           24.836474, 26.351765, 27.830157, 29.313052, 30.842023, 32.458105, 34.176166, 35.973609,
                           37.823503, 39.679749, 41.484257, 43.217056, 44.897818, 46.548041, 48.190997, 49.868754,
                           51.628845, 53.471269, 55.372669, 57.309528, 59.258842, 61.205315, 63.129999, 65.005838,
                           66.805622, 68.522718, 70.199051, 71.883548, 73.625107, 75.472584, 77.465759, 79.562769,
                           81.678562, 83.727691, 85.641583, 87.428105, 89.116967, 90.737882, 92.320415, 93.890088,
                           95.467992, 97.074803, 98.716412, 100.373750, 102.025423, 103.659693, 105.283606, 106.911009,
                           108.578746, 110.330787, 112.211101, 114.245949, 116.391751, 118.586155, 120.760791,
                           122.870526, 124.926778, 126.948832, 128.955974, 130.966890, 132.996551, 135.058676,
                           137.166984, 139.334971, 141.567627, 143.858199, 146.198901, 148.582386, 151.025877,
                           153.595036, 156.374536, 159.429517, 162.680645, 165.971583, 169.170612, 172.338348,
                           175.579778, 178.986959, 182.693415, 187.583565, 194.084702, 200.230274, 207.941894,
                           218.317735, 231.574791, 247.490165, 266.301493, 290.599491, 323.030286, 375.493924,
                           452.106557, 541.070656, 642.203415, 755.414777, 884.661680, 1025.878182, 1179.995157,
                           1350.264491, 1530.184613, 1725.812093, 1932.246478, 2151.668557, 4487.779571, 4817.359564,
                           5162.118156, 5518.688361, 5889.109108, 6276.101154, 6675.154003, 7089.259433, 7519.875831,
                           7967.994972, 8433.658528, 8921.693012, 9425.447921, 9952.111116, 10495.003204, 11066.205919,
                           11666.571692, 12286.828285, 12938.387301, 13620.164654, 14335.222658, 15079.186135,
                           15860.419436, 16682.086511, 17539.204677, 18441.998466, 19385.531299, 20375.519626,
                           21412.438530]

        self.snapshot.add_material(air)
        self.snapshot.add_material(magnet_up)
        self.snapshot.add_material(magnet_down)
        self.snapshot.add_material(printed_steel)

    def define_boundary_conditions(self):
        a0 = DirichletBoundaryCondition("a0", field_type="magnetic", magnetic_potential=0.0)

        # Adding boundary conditions to the snapshot
        self.snapshot.add_boundary_condition(a0)

    def add_postprocessing(self):
        points = [(0, 0)]
        self.snapshot.add_postprocessing("integration", points, "Energy")

    def draw_lower_line(self):
        ll = Node(1, 1)
        lr = Node(141, 1)

        O = 1 + 3.75 + 10 / 2
        w = 4

        C1 = (-0.5, 1 + 3)
        C2 = (0.5, 1 + 3)

        start_pt = []
        end_pt = []
        for i in range(8):
            start_pt.append(O - w / 2)
            end_pt.append(O + w / 2)
            c1x = C1[0] + O
            c2x = C2[0] + O

            bz = ParametricBezier(start_pt=(O - w / 2, 1),
                                  c1=(c1x, C1[1]),
                                  c2=(c2x, C2[1]),
                                  end_pt=(O + w / 2, 1))
            lx, ly = bz.approximate(4)
            for i in range(1, len(lx)):
                x0 = lx[i - 1]
                y0 = ly[i - 1]
                x1 = lx[i]
                y1 = ly[i]
                n0 = Node(lx[i - 1], ly[i - 1])
                n1 = Node(lx[i], ly[i])
                self.geom.add_line(Line(n0, n1))

            O = O + 17.5

        for i in range(len(start_pt) - 1):
            n0 = Node(start_pt[i + 1], 1)
            n1 = Node(end_pt[i], 1)
            self.geom.add_line(Line(n0, n1))

        self.geom.add_line(Line(ll, Node(start_pt[0], 1)))
        self.geom.add_line(Line(lr, Node(end_pt[-1], 1)))

    def draw_upper_line(self):
        ul = Node(6, 20)
        ur = Node(146, 20)

        O = 6 + 3.75 + 10 / 2
        w = 4

        C1 = (-0.5, 20 - 3)
        C2 = (0.5, 20 - 3)

        start_pt = []
        end_pt = []
        for i in range(8):
            start_pt.append(O - w / 2)
            end_pt.append(O + w / 2)
            c1x = C1[0] + O
            c2x = C2[0] + O

            bz = ParametricBezier(start_pt=(O - w / 2, 20),
                                  c1=(c1x, C1[1]),
                                  c2=(c2x, C2[1]),
                                  end_pt=(O + w / 2, 20))
            lx, ly = bz.approximate(4)
            for i in range(1, len(lx)):
                x0 = lx[i - 1]
                y0 = ly[i - 1]
                x1 = lx[i]
                y1 = ly[i]
                n0 = Node(lx[i - 1], ly[i - 1])
                n1 = Node(lx[i], ly[i])
                self.geom.add_line(Line(n0, n1))

            O = O + 17.5

        for i in range(len(start_pt) - 1):
            n0 = Node(start_pt[i + 1], 20)
            n1 = Node(end_pt[i], 20)
            self.geom.add_line(Line(n0, n1))

        self.geom.add_line(Line(ul, Node(start_pt[0], 20)))
        self.geom.add_line(Line(ur, Node(end_pt[-1], 20)))

    def build_geometry(self):
        g = ModelPiece('geom')
        g.load_piece_from_dxf(ModelDir.RESOURCES / 'geom.dxf')
        self.geom.merge_geometry(g.geom)

        c1l = (-1, 1)
        c2l = (1, 0.5)
        w_lower = 3
        offset_lower = (1, 1)
        
        c1u = (-1, -1)
        c2u = (1, -0.5)
        w_upper = 5
        offset_upper = (6, 20)

        for i in range(1, 9):

            origin_lower = offset_lower[0] + i * 17.5 - 7.5 # ?
            start_lower = Node(origin_lower-w_lower / 2, offset_lower[1])
            c1_lower = Node(origin_lower+c1l[0], offset_lower[1] + c1l[1])
            c2_lower = Node(origin_lower+c2l[0], offset_lower[1] + c2l[1])
            end_lower = Node(origin_lower+w_lower / 2, offset_lower[1])
            bz_lower = ParametricBezier(start_lower, c1_lower, c2_lower, end_lower)

            for li in bz_lower.approximate(n_segment=13):
                self.geom.add_line(li)

            # start_upper = Node(-w_upper / 2, offset_upper[1])
            # c1_upper = Node(c1l[0], offset_upper[1] + c1l[1])
            # c2_upper = Node(c2l[0], offset_upper[1] + c2l[1])
            # end_upper = Node(w_upper / 2, offset_upper[1])
            # bz_upper = ParametricBezier(start_upper, c1_upper, c2_upper, end_upper)

            # for li in bz_upper.approximate(n_segment=13):
            #     self.geom.add_line(li)

        self.snapshot.add_geometry(self.geom)


if __name__ == "__main__":
    m = MagneticClutch(exportname="dev")
    print(m(cleanup=False, devmode=True))
