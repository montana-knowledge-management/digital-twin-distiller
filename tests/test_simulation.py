import unittest
from digital_twin_distiller.simulationproject import SimulationProject, sim
from digital_twin_distiller.__main__ import new
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.utils import purge_dir
from pathlib import Path

CURRENT = Path(__file__).parent
MODELNAME = 'TestModel'
MODELPATH = CURRENT / MODELNAME

DEFAULT_REQ = {
    "simulation": {
        "type": "default"
    },
    "model": {},
    "tolerances": {
        "type": "ff",
        "parameters": {"x0": 0.5},
        "variables": ["T"]
    },
    "misc": {
        "processes": 4,
        "cleanup": True
    },
    "version": "0.7"
}


class TestSimulation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # creates a dummy model
        new(MODELNAME, CURRENT)

        # import the model class from the new model
        cls.modelclass = __import__(str('tests.' + MODELNAME + '.model'), fromlist=['TestModel']).TestModel

        # set the paths for the new model
        ModelDir.set_base(MODELPATH)

    @classmethod
    def tearDownClass(cls):
        pass
        # CLEANUP SECTION

        # clean up the modeldir
        # purge_dir(MODELPATH) # DO NOT MODIFY THIS LINE

    def test_init(self):
        sim1 = SimulationProject(None)
        self.assertEqual(sim1.model, None)

    def test_set_model(self):
        sim1 = SimulationProject(None)
        sim1.set_model(self.modelclass)
        self.assertEqual(type(self.modelclass), type(self.modelclass))

    def test_load_defaults(self):
        sim1 = SimulationProject(self.modelclass)
        sim1._input = DEFAULT_REQ
        sim1._load_defaults()

        self.assertAlmostEqual(sim1.cfg_model['x0'], 1.0, delta=1e-12)
        self.assertAlmostEqual(sim1.cfg_model['mw'], 5.0, delta=1e-12)

        # Default simulation parameters
        self.assertAlmostEqual(sim1.cfg_simulation['t0'], 0.0, delta=1e-12)
        self.assertAlmostEqual(sim1.cfg_simulation['t1'], 5.3, delta=1e-12)
        self.assertEqual(int(sim1.cfg_simulation['nstep']), 101)

        self.assertFalse(sim.cfg_tolerances)

        self.assertEqual(sim1.cfg_misc['processes'], 4)
        self.assertEqual(sim1.cfg_misc['cleanup'], True)

    def test_update_input(self):
        sim1 = SimulationProject(self.modelclass)
        sim1._input = DEFAULT_REQ
        sim1._load_defaults()

        sim1.update_input()

        self.assertAlmostEqual(sim1.cfg_tolerances["parameters"]["x0"], 0.5, delta=1e-12)
        self.assertTrue('T' in sim1.cfg_tolerances['variables'])
        self.assertTrue(sim1.cfg_tolerances['type'], 'ff')


    # def test_run(self):
    #     sim1 = SimulationProject(self.modelclass)
    #     sim1._input = DEFAULT_REQ
    #     sim1._load_defaults()
    #
    #     sim1.update_input()
