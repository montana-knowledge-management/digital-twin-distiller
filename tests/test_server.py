import unittest
from multiprocessing import Process
from pathlib import Path
import os

from digital_twin_distiller.__main__ import new
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.server import Server
from digital_twin_distiller.simulationproject import sim, SimulationProject
from digital_twin_distiller.utils import purge_dir
from time import sleep
import requests

CURRENT = Path(__file__).parent
MODELNAME = "TestModel"
MODELPATH = CURRENT / MODELNAME

# class TestIntegratedServer(unittest.TestCase):
#     serverprocess = None
#
#     @classmethod
#     def setUpClass(cls):
#         # creates a dummy model
#         new(MODELNAME, CURRENT)
#
#         # import the model class from the new model
#         modelclass = __import__(str("tests." + MODELNAME + ".model"), fromlist=["TestModel"]).TestModel
#
#         # set the paths for the new model
#         ModelDir.set_base(MODELPATH)
#
#         # plug the model into the simulation
#         sim.set_model(modelclass)
#
#         # plug the simulation into the server
#         model = Server(sim)
#
#         # Create a new Process object
#         cls.serverprocess = Process(name="testserver_process", target=model, daemon=True)
#
#         # Fire up the server in a new process
#         cls.serverprocess.start()
#
#     @classmethod
#     def tearDownClass(cls):
#         # CLEANUP SECTION
#
#         # kill the serverprocess
#         cls.serverprocess.kill()
#
#         # clean up the modeldir
#         purge_dir(MODELPATH)  # DO NOT MODIFY THIS LINE TODO: OT: why???
#
#     # def test_ping(self):
#     #     sleep(2)
#     #     url = "http://0.0.0.0:5000"
#     #     res = requests.get(f'{url}/ping', timeout=10, )
#     #     print(res.json())
#     #     self.assertEqual(res.status_code, 200)
#
#     def test_root(self):
#         sleep(2)
#         url = "http://0.0.0.0:5000"
#         res = requests.get(f'{url}/', timeout=10, )
#         print(res.json())
#         self.assertEqual(res.status_code, 200)
#
# # def test_missing_text_key_sim(self):
# #     sleep(2)
# #     wrong_json = {"test": "Without Text key."}
# #     url = "http://0.0.0.0:5000"
# #     response = requests.post(f'{url}/process_sim', timeout=10, json=wrong_json,
# #                              headers={"Content-Type": "application/json"})
# #     print(response.json())
# #     self.assertTrue(response.status_code, 200)
# #     self.assertDictEqual(
# #         response.json(),
# #         {"detail": [{"loc": ["body", "simulation"], "msg": "field required", "type": "value_error.missing"}]},
# #     )
#
# # TESTING SECTION
# # Testing goes here, use requests library to pass requests to the server.
#
# # TODO: GK: Temporary placeholder
# # self.assertTrue(True)
#
# #
# # from digital_twin_distiller.server import Server
# #
# # model = Server(sim)
# # model.build_docs()
# # model.run()
#
# # for path in Path("test_server").glob("**/*"):
# #     if path.is_file():
# #         path.unlink()
# #     elif path.is_dir():
# #         rmtree(path)
# # rmtree('test_server')


from digital_twin_distiller.ml_project import MachineLearningProject
from fastapi.testclient import TestClient


class DummyMLandSimulationProject(MachineLearningProject):
    _input = []
    _output = []

    def run(self):
        self._output_data = self._input_data
        self._output = self._input

    def update_input(self):
        pass


class TestServer(unittest.TestCase):
    new(MODELNAME, CURRENT)
    # HTTP client
    example_project = DummyMLandSimulationProject(app_name="test_name")
    server = Server(example_project)
    # server.set_project_mkdocs_dir_path(ModelDir.DOCS)
    client = TestClient(server.app)

    # HTTPS client
    server2 = Server(example_project)
    server2.set_key_file_path("some_path")
    server2.set_cert_file_path("some_path")

    # https_client = TestClient(server2.app)

    def test_http_ping(self):
        response = self.client.get("/ping")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["msg"], "The API is working.")
        self.assertIn("call_time", response.json())

    # def test_http_root(self):
    #     response = self.client.get("/")
    #     self.assertEqual(response.status_code, 200)

    def test_server(self):
        self.assertEqual(self.server.host, "127.0.0.1")
        self.assertEqual(self.server.port, 5000)
        self.assertEqual(self.server.cert_file_path, None)
        self.assertEqual(self.server.key_file_path, None)
        self.assertEqual(self.server.app.title, "test_name API")

    def test_missing_text_key_ml(self):
        wrong_json = {"test": "Without Text key."}
        response = self.client.post("/process_ml", json=wrong_json, headers={"Content-Type": "application/json"})
        self.assertIsNot(response.status_code, 200)
        self.assertDictEqual(
            response.json(),
            {"detail": [{"loc": ["body", "text"], "msg": "field required", "type": "value_error.missing"}]},
        )

    def test_with_text_key_sim(self):
        response = self.client.post("/process_sim", json={}, headers={"Content-Type": "application/json"})
        self.assertTrue(response.status_code, 200)
        self.assertDictEqual(response.json(), {'simulation': {'type': 'default'}, 'model': {},
                                               'tolerances': {'type': 'ff', 'parameters': {}, 'variables': []},
                                               'misc': {'processes': 4, 'cleanup': True, 'exportname': None},
                                               'version': '2021.12'}
                             )

    def test_with_text_key_ml(self):
        good_json = {"text": "With text key."}
        response = self.client.post("/process_ml", json=good_json, headers={"Content-Type": "application/json"})
        self.assertTrue(response.status_code, 200)
        self.assertDictEqual(response.json(), good_json)

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

    def test_asset_not_exists(self):
        with self.assertRaises(FileNotFoundError) as context:
            self.server.set_project_mkdocs_dir_path(".")
            self.assertIn('please build the mkdocs site by calling "mkdocs build" in the docs directory',
                          context.exception)
