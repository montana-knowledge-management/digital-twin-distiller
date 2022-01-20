from unittest import TestCase

from importlib_resources import files

from digital_twin_distiller.geo_parser import geo_parser


class SimpleGEOparser(TestCase):
    def test_base_objects(self):
        path = files("tests.test_geo_geometry").joinpath("base_objects.geo")
        test_geo = geo_parser(path.as_posix())

        # the point data and the coordinates in the given geo file
        self.assertEqual(test_geo.nodes[0].x, -0.0004018925337573828)
        self.assertEqual(test_geo.nodes[0].y, 0.05177306826503763)
        self.assertEqual(test_geo.nodes[0].id, 1336)

        # the line object imported correctly
        self.assertEqual(test_geo.lines[0].id, 1)
        self.assertEqual(test_geo.lines[0].start_pt.id, 1)
        self.assertEqual(test_geo.lines[0].end_pt.id, 2)

        # import a circle object circle data 60 [1.0, 2.0, 1336.0]
        self.assertEqual(test_geo.circle_arcs[0].id, 60)
        self.assertEqual(test_geo.circle_arcs[0].start_pt.id, 1336)
        self.assertEqual(test_geo.circle_arcs[0].end_pt.id, 1)

        # import second circle 160 [101, 102, 103]
        self.assertEqual(test_geo.circle_arcs[1].id, 160)
        self.assertEqual(test_geo.circle_arcs[1].start_pt.id, 103)
        self.assertEqual(test_geo.circle_arcs[1].center_pt.id, 102)
        self.assertEqual(test_geo.circle_arcs[1].end_pt.id, 101)

        # point 3 should be imported from the geo file
        self.assertEqual(test_geo.nodes[3].id, 3)
        self.assertEqual(test_geo.nodes[3].x, 0.1)

    def test_base_object_invalid(self):
        path = files("tests.test_geo_geometry").joinpath("invalid_line_object.geo")
        test_geo = geo_parser(path.as_posix())
        # this test should handle the comments properly and shouldnt fail on the example file
        # the only one valid line object should  imported correctly
        self.assertEqual(test_geo.lines[0].id, 1)
        self.assertEqual(len(test_geo.lines), 1)
        self.assertEqual(test_geo.lines[0].end_pt.id, 2)

    def test_first_gmsh_geo_example(self):
        path = files("tests.test_geo_geometry").joinpath("geo_import_test.geo")
        test_geo = geo_parser(path.as_posix())

        self.assertEqual(test_geo.lines[0].id, 1)
        self.assertEqual(len(test_geo.lines), 4)
