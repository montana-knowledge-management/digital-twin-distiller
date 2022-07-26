import unittest
from pathlib import Path

import requests
from fastapi.testclient import TestClient

from digital_twin_distiller.encapsulator import Encapsulator
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import SimulationProject
from digital_twin_distiller.utils import purge_dir
from digital_twin_distiller.cli import new

CURRENT = Path(__file__).parent
MODELNAME = "TestServer"
MODELPATH = CURRENT / MODELNAME


class DummySimulationProject(SimulationProject):
    _input = []
    _output = []

    def run(self):
        self._output = self._input

    def update_input(self):
        pass


class TestServer(unittest.TestCase):
    new(MODELNAME, CURRENT)

    # HTTP client
    example_project = DummySimulationProject()
    server = Encapsulator(example_project)
    server.set_endpoint_for_docs(ModelDir.DOCS)
    client = TestClient(server.app)

    # HTTPS client
    server2 = Encapsulator(example_project)
    server2.set_key_file_path("some_path")
    server2.set_cert_file_path("some_path")

    def test_http_ping(self):
        response = self.client.get("/ping")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["msg"], "The API is working.")
        self.assertIn("call_time", response.json())

    def test_http_root(self):
        response = self.client.get("/")
        self.assertEqual(response.status_code, 200)

    def test_server(self):
        self.assertEqual(self.server.host, "127.0.0.1")
        self.assertEqual(self.server.port, 5000)
        self.assertEqual(self.server.cert_file_path, None)
        self.assertEqual(self.server.key_file_path, None)
        self.assertEqual(self.server.app.title, "digital twin project API")

    def test_with_empty_input_sim(self):
        response = self.client.post("/process_sim", json={}, headers={"Content-Type": "application/json"})
        self.assertTrue(response.status_code, 200)
        self.assertDictEqual(
            response.json(),
            {
                "simulation": {"type": "default"},
                "model": {},
                "tolerances": {"type": "ff", "parameters": {}, "variables": []},
                "misc": {"processes": 4, "cleanup": True, "exportname": None},
                "version": "2021.12",
            },
        )

    def test_https_paths(self):
        self.assertIsNotNone(self.server2.cert_file_path)
        self.assertIsNotNone(self.server2.key_file_path)

    def test_set_host(self):
        example_host = "1.2.3.4"
        # serv = Server(self.example_project)
        self.server.set_host(example_host)
        self.assertEqual(self.server.host, example_host)

    def test_set_port(self):
        example_port = 123
        self.server.set_port(example_port)
        self.assertEqual(self.server.port, example_port)

    @classmethod
    def tearDownClass(cls):
        # CLEANUP SECTION
        # clean up the modeldir
        purge_dir(MODELPATH)
