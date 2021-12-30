from unittest import TestCase
from examples.fem_simulations.power_transformer.model import PowerTransformer


class TestPowerTransformerExample(TestCase):
    """
    The tests performed with the data of an existing transformer: DOI: 10.3233/JAE-209504
    --------------------------------------------------------------------------------------
    Nominal power MVA 10
    Frequency Hz 50
    Connection group YNyn1
    Number of phases # 3
    Short circuit impedance % 7.5
    Main gap mm 50
    Sum of the end insulation mm 80
    Phase distance mm 50
    Core-Inner winding distance mm 20
    Number of legs # 3
    Flux density limit in columns T 1.7
    Core Filling Factor % 88
    Material Type M1H
    Material Price e/kg 3.5
    Line Voltage kV 6.9
    Phase Voltage kV 4.0
    Low Voltage BIL kV 125
    Winding AC kV 50
    Copper filling factor % 60
    Material and manufacturing price e/kg 10
    Line Voltage kV 33
    Phase Voltage kV 19.05
    High Voltage BIL kV 125
    Winding AC kV 50
    Copper filling factor % 60
    Material and manufacturing price e/kg 8.5
    """

    def init_transformer_model(self):
        design_parameters = {'ff_in': 60,  # filling factor in the inner winding [%]
                             'ff_ou': 60,  # filling factor in the outer winding [%]
                             'alpha': 1.0,  # ratio
                             'end_ins': 80,  # the sum of the end insulation under and on the windings [mm]
                             'core_ins': 20}  # core insulation distance between the core and the windings [mm]

        design_variables = {'core_diam': 420,  # core diameter in mm
                            'gap': 50,  # main insulation distance between the main windings in mm
                            'hin': 1100,  # the height of the inner winding
                            'tin': 35,  # thickness of the inner winding
                            'tou': 44,  # thickness of the outer winding
                            'jin': 2.6,  # the current density in the inner winding
                            'jou': 2.5}  # the current density in the outer winding

        self.pt = PowerTransformer(params=design_parameters, vars=design_variables)

    def test_init_transformer_model(self):
        self.init_transformer_model()

        self.assertEqual(self.pt.design_parameters['ff_in'], 60)
        self.assertEqual(self.pt.design_parameters['ff_ou'], 60)
        self.assertEqual(self.pt.design_parameters['end_ins'], 80)
        self.assertEqual(self.pt.design_parameters['alpha'], 1.0)
        self.assertEqual(self.pt.design_parameters['core_ins'], 20)

        self.assertEqual(self.pt.design_variables['core_diam'], 420)
        self.assertEqual(self.pt.design_variables['gap'], 50)
        self.assertEqual(self.pt.design_variables['jin'], 2.6)
        self.assertEqual(self.pt.design_variables['jou'], 2.5)

    def test_calculate_window_parameters(self):
        self.init_transformer_model()
        self.pt.set_window_parameters()

        self.assertEqual(self.pt.w1, 199)
        self.assertEqual(self.pt.h1, 1180)

        #
        self.assertEqual(self.pt.r1, 210)
        self.assertEqual(self.pt.z1, 40)

    def test_windings(self):
        self.init_transformer_model()
        self.pt.set_window_parameters()
        self.pt.set_inner_winding_parameters()
        self.pt.set_outer_winding_parameters()

        self.assertEqual(self.pt.w2, 35)
        self.assertEqual(self.pt.h2, 1100)
        self.assertEqual(self.pt.r2, 230)
        self.assertEqual(self.pt.z2, 40)

        self.assertEqual(self.pt.w3, 44)
        self.assertEqual(self.pt.h3, 1100)
        self.assertEqual(self.pt.r3, 315)
        self.assertEqual(self.pt.z3, 40)
