from unittest import TestCase
from examples.fem_simulations.power_transformer.model import PowerTransformer
from examples.fem_simulations.power_transformer.simulation import calculate_base_impedance
from math import pi


class TestPowerTransformerExample(TestCase):
    """
    The tests performed with the data of an existing transformer: DOI: 10.3233/JAE-209504
    --------------------------------------------------------------------------------------
    Turns in LV winding 85
    area of one conductor=22.684×10^-6 m2
    """

    def init_transformer_model(self):
        self.pt = PowerTransformer(ff_in=60, ff_ou=60, alpha=1.0, end_ins=180., core_ins=20., core_diam=420., gap=50.,
                                   hin=1100, tin=35, tou=44, jin=2.6, jou=2.26, exportname="dev")

    def test_init_transformer_model(self):
        self.init_transformer_model()

        self.assertEqual(self.pt.design_parameters['ff_in'], 60)
        self.assertEqual(self.pt.design_parameters['ff_ou'], 60)
        self.assertEqual(self.pt.design_parameters['end_ins'], 180)
        self.assertEqual(self.pt.design_parameters['alpha'], 1.0)
        self.assertEqual(self.pt.design_parameters['core_ins'], 20)

        self.assertEqual(self.pt.design_variables['core_diam'], 420)
        self.assertEqual(self.pt.design_variables['gap'], 50)
        self.assertEqual(self.pt.design_variables['jin'], 2.6)
        self.assertEqual(self.pt.design_variables['jou'], 2.26)

    def test_calculate_window_parameters(self):
        self.init_transformer_model()

        self.assertEqual(self.pt.w1, 199)
        self.assertEqual(self.pt.h1, 1280)

        self.assertEqual(self.pt.r1, 210)
        self.assertEqual(self.pt.z1, 0)

    def test_windings(self):
        self.init_transformer_model()

        self.assertEqual(self.pt.w2, 35)
        self.assertEqual(self.pt.h2, 1100)
        self.assertEqual(self.pt.r2, 230)
        self.assertEqual(self.pt.z2, 90)
        self.assertEqual(self.pt.js, 1560000.0)

        self.assertEqual(self.pt.w3, 44)
        self.assertEqual(self.pt.h3, 1100)
        self.assertEqual(self.pt.r3, 315)
        self.assertEqual(self.pt.z3, 90)
    #  TODO
    # def test_total_model(self):
    #     self.init_transformer_model()
    #
    #     # calculating the magnetic energy
    #     Wm = self.pt(cleanup=False, devmode=False)['Energy']
    #
    #     # base impedance for the short circuit impedance calculation
    #     ub = 6.9  # voltage --- kV base voltage
    #     sb = 10.0  # nominal power  --- MVA
    #     f = 50
    #
    #     sci = calculate_base_impedance(ub, sb, f, 3 ** 0.5, Wm)
    #
    #     print('Short circuit impedance', sci)
    #     self.assertEqual(7.4, round(sci, 1))
