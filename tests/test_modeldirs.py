import unittest
from digital_twin_distiller.modelpaths import ModelDir
from pathlib import Path

CURRENT = Path(__file__)

class TestModeldirs(unittest.TestCase):

    def test_modeldir(self):
        ModelDir.set_base(CURRENT)

        self.assertEqual(ModelDir.BASE, CURRENT.parent)
        self.assertEqual(ModelDir.MEDIA, CURRENT.parent / "media")
        self.assertEqual(ModelDir.DATA, CURRENT.parent / "data")
        self.assertEqual(ModelDir.RESOURCES, CURRENT.parent / "resources")
        self.assertEqual(ModelDir.SNAPSHOTS, CURRENT.parent / "snapshots")
        self.assertEqual(ModelDir.DEFAULTS, CURRENT.parent / "defaults")
        self.assertEqual(ModelDir.DOCS, CURRENT.parent / "docs")


    def test_dir_iter(self):
        ModelDir.set_base(CURRENT.parent)

        dirs = ModelDir.get_dirs()

        self.assertEqual(next(dirs), CURRENT.parent)
        self.assertEqual(next(dirs), CURRENT.parent / "media")
        self.assertEqual(next(dirs), CURRENT.parent / "data")
        self.assertEqual(next(dirs), CURRENT.parent / "resources")
        self.assertEqual(next(dirs), CURRENT.parent / "snapshots")
        self.assertEqual(next(dirs), CURRENT.parent / "defaults")
        self.assertEqual(next(dirs), CURRENT.parent / "docs")
