import argparse
import io
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

from digital_twin_distiller import cli
from digital_twin_distiller import purge_dir
from digital_twin_distiller.cli import NAME_OF_THE_PROGRAM
from digital_twin_distiller.cli import optimize_cli

CURRENT = Path(__file__).parent

TEST_NAME = "cli_test_name"
TEST_FOLDER = "cli_test_folder"
MODEL_DIR = CURRENT / TEST_FOLDER


class TestCli(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        cli.NAME_OF_THE_PROGRAM = NAME_OF_THE_PROGRAM

    @patch('argparse.ArgumentParser.parse_args',
           return_value=argparse.Namespace(command="new", name=TEST_NAME, location=MODEL_DIR))
    def test_parse_arguments(self, mock_args):
        """testing the required params (command, new, location)"""
        optimize_cli(["new", TEST_NAME, MODEL_DIR])
        purge_dir(MODEL_DIR, force=True)

    def test_new_model_creation(self):
        """testing if the model template created by calling 'new' representing from command line"""
        optimize_cli(["new", TEST_NAME, str(MODEL_DIR)])
        self.assertTrue(MODEL_DIR.exists())
        purge_dir(MODEL_DIR, force=True)

    def test_invalid_param_list_when_new_called(self):
        stderr = io.StringIO()
        sys.stderr = stderr
        with self.assertRaises(SystemExit):
            optimize_cli(["new"])
        self.assertRegexpMatches(stderr.getvalue(), r"the following arguments are required")

    @unittest.skip("output not working after 'test_invalid_param_list_when_new_called' has been called")
    def test_when_required_param_is_unknown(self):
        """ Try to perform when param isn't an option. """
        stderr = io.StringIO()
        sys.stderr = stderr
        with self.assertRaises(SystemExit):
            optimize_cli(["unknown"])
        self.assertRegexpMatches(stderr.getvalue(), r"invalid choice")

    def test_output_when_version_called(self):
        self.valid_version_call(param="-v")
        self.valid_version_call(param="--version")

    def valid_version_call(self, param: str):
        stored_out = io.StringIO()
        sys.stdout = stored_out
        with self.assertRaises(SystemExit):
            optimize_cli([param])
        self.assertRegexpMatches(stored_out.getvalue(), r"digital-twin-distiller [\d.]+ \nPython [\d.]+")

    def test_when_exception_thrown(self):
        cli.NAME_OF_THE_PROGRAM = "unknown-package-name"
        with self.assertRaises(SystemExit) as cm:
            optimize_cli()
        self.assertEqual(1, cm.exception.code)
