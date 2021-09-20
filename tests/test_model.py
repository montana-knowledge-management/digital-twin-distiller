import unittest
from adze_modeler.metadata import Agros2DMetadata
from adze_modeler.model import BaseModel
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.snapshot import Snapshot
from pathlib import Path
from shutil import rmtree


class MockModel(BaseModel):
    def __init__(self, exportname: str = None):
        super().__init__(exportname)

    def setup_solver(self):
        super().setup_solver()
        # Agros2D
        agros_metadata = Agros2DMetadata()
        agros_metadata.file_script_name = self.file_solver_script
        agros_metadata.file_metrics_name = self.file_solution
        agros_metadata.problem_type = "magnetic"
        agros_metadata.coordinate_type = "axisymmetric"
        agros_metadata.analysis_type = "steadystate"
        agros_metadata.unit = 1e-3
        agros_metadata.nb_refinements = 0
        agros_metadata.adaptivity = "hp-adaptivity"
        agros_metadata.adaptivity_tol = 1
        platform_agros = Agros2D(agros_metadata)
        self.snapshot = Snapshot(platform_agros)

    def add_postprocessing(self):
        super().add_postprocessing()

    def define_materials(self):
        super().define_materials()

    def define_boundary_conditions(self):
        super().define_boundary_conditions()

    def build_geometry(self):
        super().build_geometry()


class TestModel(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dir_current = Path(__file__).parent
        self.dir_resources = self.dir_current / "resources"
        self.dir_snapshots = self.dir_current / "snapshots"
        self.dir_media = self.dir_current / "media"
        self.dir_data = self.dir_current / "data"
        # self.dir_export = self.dir_snapshots / self.name
        self.dir_export = self.dir_snapshots

    def test_model_paths(self):
        m = MockModel(exportname="test_name")
        self.assertEqual(self.dir_current, m.dir_current)
        self.assertEqual(self.dir_resources, m.dir_resources)
        self.assertEqual(self.dir_snapshots, m.dir_snapshots)
        self.assertEqual(self.dir_media, m.dir_media)
        self.assertEqual(self.dir_data, m.dir_data)
        self.assertEqual(self.dir_export, m.dir_export.parent)

        # create paths
        m._init_directories()
        self.assertTrue(self.dir_current.exists())
        self.assertTrue(self.dir_resources.exists())
        self.assertTrue(self.dir_snapshots.exists())
        self.assertTrue(self.dir_media.exists())
        self.assertTrue(self.dir_data.exists())
        self.assertTrue(self.dir_export.exists())

        m(cleanup=True, timeout=0)
        rmtree(m.dir_snapshots)

        # teardown
        self.dir_resources.rmdir()
        self.dir_media.rmdir()
        self.dir_data.rmdir()
