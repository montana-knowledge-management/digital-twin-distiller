import argparse
import unittest
from unittest.mock import patch

import digital_twin_distiller
from digital_twin_distiller.cli import optimize_cli
from unittest import mock
from digital_twin_distiller.__main__ import new


class TestCli(unittest.TestCase):

    @mock.patch('argparse.ArgumentParser.parse_args',
                return_value=argparse.Namespace(command="new", name="test_name", location="sandbox"))
    def test_parse_arguments(self, mock_args):
        optimize_cli(["new", "test_name", "sandbox"])

    @patch.object(digital_twin_distiller.__main__, "new")
    def test_when_new_method_is_called(self, mock_args):
        # digital_twin_distiller.cli.optimize_cli(["new", "test_name", "sandbox"])
        t = digital_twin_distiller.cli
        t.new("test_name", "sandbox")
        t.optimize_cli(["new", "a", "b"])
        unittest.mock.Mock.assert_called(mock_args)



    def test_greet_cli(capsys):
        optimize_cli(["-h"])

        captured = capsys.readouterr()

        result = captured.out
        assert result.find("Welcome to Digital Twin Distiller!")
