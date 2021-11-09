import unittest
from pathlib import Path
import requests
from time import sleep

from digital_twin_distiller.__main__ import new
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import sim
from digital_twin_distiller.server import Server
from multiprocessing import Process
from digital_twin_distiller.utils import purge_dir

CURRENT = Path(__file__).parent
MODELNAME = 'TestModel'
MODELPATH = CURRENT / MODELNAME

class TestIntegratedServer(unittest.TestCase):
    serverprocess=None

    @classmethod
    def setUpClass(cls):
        # creates a dummy model
        new(MODELNAME, CURRENT)

        # import the model class from the new model
        modelclass = __import__(str('tests.'+MODELNAME+'.model'), fromlist=['TestModel']).TestModel

        # set the paths for the new model
        ModelDir.set_base(MODELPATH)

        # plug the model into the simulation
        sim.set_model(modelclass)

        # plug the simulation into the server
        model = Server(sim)

        # Create a new Process object
        cls.serverprocess = Process(name='testserver_process',
                                    target=model,
                                    daemon=True)

        # Fire up the server in a new process
        cls.serverprocess.start()

    @classmethod
    def tearDownClass(cls):
        # CLEANUP SECTION

        # kill the serverprocess
        cls.serverprocess.kill()

        # clean up the modeldir
        purge_dir(MODELPATH) # DO NOT MODIFY THIS LINE

    def test_ping(self):
        # sleep(3)
        url = "http://0.0.0.0:5000"
        res = requests.get(f'{url}/apidocs', timeout=10,)
        self.assertEqual(res.status_code, 200)


        # TESTING SECTION
        # Testing goes here, use requests library to pass requests to the server.
























#         #
#         # from digital_twin_distiller.server import Server
#         #
#         # model = Server(sim)
#         # model.build_docs()
#         # model.run()
#
#         # for path in Path("test_server").glob("**/*"):
#         #     if path.is_file():
#         #         path.unlink()
#         #     elif path.is_dir():
#         #         rmtree(path)
#         rmtree('test_server')


# from digital_twin_distiller.server import Server
# from examples.text_classification.deployment.cached_twenty_news_project import CachedExampleProject
# from fastapi.testclient import TestClient
#
# # HTTP client
# example_project = CachedExampleProject(app_name="test_name")
# server = Server(example_project)
# client = TestClient(server.app)
#
# # HTTPS client
# server2 = Server(example_project)
# server2.set_key_file_path("some_path")
# server2.set_cert_file_path("some_path")
# https_client = TestClient(server2.app)
#
#
# class TestIntegratedServer(unittest.TestCase):
#     def test_http_ping(self):
#         response = client.get("/ping")
#         self.assertEqual(response.status_code, 200)
#         self.assertEqual(response.json()["msg"], "The API is working.")
#         self.assertIn("call_time", response.json())
#
#     def test_http_root(self):
#         response = client.get("/")
#         self.assertEqual(response.status_code, 200)
#
#     def test_server(self):
#         self.assertEqual(server.host, "127.0.0.1")
#         self.assertEqual(server.port, 5000)
#         self.assertEqual(server.cert_file_path, None)
#         self.assertEqual(server.key_file_path, None)
#         self.assertEqual(server.app.title, "test_name API")
#
#     def test_http_post(self):
#         input_json = {"text": "This is a Christian test.", "ExtraKey": "Extra."}
#         response = client.post("/process", json=input_json, headers={"Content-Type": "application/json"})
#         self.assertEqual(response.status_code, 200)
#         self.assertEqual(response.json().get("Result"), "soc.religion.christian")
#         self.assertEqual(input_json.get("text"), response.json().get("text"))
#         self.assertEqual(input_json.get("ExtraKey"), response.json().get("ExtraKey"))
#
#     def test_missing_text_key(self):
#         wrong_json = {"test": "Without Text key."}
#         response = client.post("/process", json=wrong_json, headers={"Content-Type": "application/json"})
#         self.assertIsNot(response.status_code, 200)
#         self.assertDictEqual(
#             response.json(),
#             {"detail": [{"loc": ["body", "text"], "msg": "field required", "type": "value_error.missing"}]},
#         )
#
#     def test_https_paths(self):
#         self.assertIsNotNone(server2.cert_file_path)
#         self.assertIsNotNone(server2.key_file_path)
