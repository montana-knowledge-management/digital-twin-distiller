import unittest
from adze_modeler.__main__ import new
from pathlib import Path
from adze_modeler.utils import purge_dir

CURRENT = Path(__file__).parent

class TestMain(unittest.TestCase):
    def test_new_model(self):
        MODELNAME = 'TestModel'
        MODELPATH = CURRENT / MODELNAME
        new(MODELNAME, CURRENT)
        mustexists = lambda item: self.assertTrue((MODELPATH / item).exists())

        mustexists('data')
        mustexists('defaults')
        mustexists('defaults/misc.json')
        mustexists('defaults/model.json')
        mustexists('defaults/simulation.json')
        mustexists('docs')
        mustexists('docs/mkdocs.yml')
        mustexists('docs/docs')
        mustexists('docs/site')
        mustexists('media')
        mustexists('resources')
        mustexists('snapshots')
        mustexists('model.py')
        mustexists('simulation.py')

        purge_dir(MODELPATH, force=True)
