import unittest
from pathlib import Path

from digital_twin_distiller.cli import new
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import SimulationProject
from digital_twin_distiller.utils import purge_dir

CURRENT = Path(__file__).parent
MODELNAME = "TestModel"
MODELPATH = CURRENT / MODELNAME

DEFAULT_REQ = {
    "simulation": {"type": "default"},
    "model": {},
    "tolerances": {"type": "ff", "parameters": {"x0": 1, "mw": 1}, "variables": ["T"]},
    "misc": {"processes": 4, "cleanup": True, "exportname": "testname"},
    "version": "0.7",
}


class TestSimulation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # creates a dummy model
        new(MODELNAME, CURRENT)

        # import the model class from the new model
        cls.modelclass = __import__(str("tests." + MODELNAME + ".model"), fromlist=[MODELNAME]).SimulationModel

        # set the paths for the new model
        ModelDir.set_base(MODELPATH)

    @classmethod
    def tearDownClass(cls):
        # CLEANUP SECTION

        # clean up the modeldir
        purge_dir(MODELPATH)  # DO NOT MODIFY THIS LINE

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

        self.assertAlmostEqual(sim1.cfg_model["x0"], 1.0, delta=1e-12)
        self.assertAlmostEqual(sim1.cfg_model["mw"], 5.0, delta=1e-12)

        # Default simulation parameters
        self.assertAlmostEqual(sim1.cfg_simulation["t0"], 0.0, delta=1e-12)
        self.assertAlmostEqual(sim1.cfg_simulation["t1"], 5.3, delta=1e-12)
        self.assertEqual(int(sim1.cfg_simulation["nstep"]), 101)

        self.assertFalse(sim1.cfg_tolerances)

        self.assertEqual(sim1.cfg_misc["processes"], 4)
        self.assertEqual(sim1.cfg_misc["cleanup"], True)

    def test_update_input(self):
        sim1 = SimulationProject(self.modelclass)
        sim1._input = DEFAULT_REQ
        sim1._load_defaults()

        sim1.update_input()

        self.assertEqual(sim1.cfg_model["exportname"], "testname")
        self.assertAlmostEqual(sim1.cfg_tolerances["parameters"]["x0"], 1.0, delta=1e-12)
        self.assertTrue("T" in sim1.cfg_tolerances["variables"])
        self.assertTrue(sim1.cfg_tolerances["type"], "ff")

    def test_run(self):
        def default_sim(model, modelparams, simparams, miscparams):
            x0 = modelparams["x0"]
            mw = modelparams["mw"]
            return {"T": 1 + x0 * 10 + mw * 100}

        sim1 = SimulationProject(self.modelclass)
        sim1._input = DEFAULT_REQ
        sim1._load_defaults()
        sim1.update_input()

        sim1.register("default")(default_sim)

        sim1.run()
        self.assertAlmostEqual(sim1._output["res"]["T"], 511.0)
