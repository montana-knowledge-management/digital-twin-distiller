import socket
import unittest

from importlib_resources import files

import digital_twin_distiller.ml_project as mlp
import digital_twin_distiller.text_readers as r
from digital_twin_distiller.keywords import JSON, PDF, TXT


class DummySubTask(mlp.AbstractSubTask):
    def run(self, input_data):
        return input_data.get("test")


class TestMLProject(unittest.TestCase):
    def test_dummy_project_ip_helper(self):
        # tests the basic functions of an empty abstractproject
        Project = mlp.MachineLearningProject()

        # init
        self.assertIsInstance(Project, mlp.MachineLearningProject)

        # check ip address helper
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
        s.connect(("<broadcast>", 0))
        ip = s.getsockname()[0]

        self.assertEqual(Project.get_ip()[0], ip)

    def test_dummy_project_text_file_reader(self):
        # tests the basic functions of an empty abstractproject
        Project = mlp.MachineLearningProject()

        # check the input_processing helpers on a text and pdf files
        text_file = files("tests") / "test_documents" / "Moore.txt"
        pdf_file = files("tests") / "test_documents" / "Moore.pdf"

        # input under a Text key
        inp = Project.open_data_file(r.TextReader(), text_file, key="Text")
        self.assertIn("Moore", inp["Text"])

        # input under a PdfText key
        inp2 = Project.open_data_file(r.PdfReader(), pdf_file, key="PdfText")
        # dictionary under a key
        self.assertIn("Moore", inp2["PdfText"]["Text"])

        inp3 = Project.open_data_file(r.PdfReader(), pdf_file)
        self.assertIn("Moore", inp3["Text"])

    def test_supported_extensions(self):
        self.assertListEqual(mlp.supported_extensions, [JSON, TXT, PDF])

    def test_input_output(self):
        example_input = {"test": "example"}
        project = mlp.MachineLearningProject()
        project.add_single_input(example_input)
        self.assertDictEqual(project._input_data[0], example_input)
        self.assertEqual(project.get_single_output(), None)
        project._output_data = ["test_data"]
        self.assertEqual(project.get_single_output(), "test_data")

    def test_abstract_task_init(self):
        example_data = "example_string"
        abstract_task = mlp.AbstractTask(project_inp=example_data)
        self.assertEqual(abstract_task._input_data, example_data)
        self.assertEqual(abstract_task._output_data, {})
        self.assertEqual(abstract_task._temp_data, {})
        self.assertEqual(abstract_task.sub_tasks, [])

    def test_simple_execution_abstract_task(self):
        onetwothree = "123"
        example_dict = {"test": onetwothree}
        abstract_task = mlp.AbstractTask(project_inp=[example_dict, example_dict])
        abstract_task.sub_tasks = [DummySubTask()]
        ret_list = abstract_task.simple_execution()
        self.assertEqual(ret_list, [onetwothree, onetwothree])
        abstract_task._input_data = [example_dict]
        result = abstract_task([DummySubTask()])
        self.assertEqual(result, [onetwothree])

    def test_subtask(self):
        data = {"test": 2}
        dummy = DummySubTask()
        self.assertEqual(dummy(data), 2)

    def test_bulk_input_directory(self):
        test_doc_path = files("tests") / "test_documents"
        project = mlp.MachineLearningProject()
        # testing pdf
        project.bulk_input_directory(directory_name=test_doc_path, extension=".pdf", key="PDF")
        self.assertIsNotNone(project._input_data[0].get("PDF"))
        self.assertEqual(len(project._input_data), 2)
        # testing txt
        project._input_data = []
        project.bulk_input_directory(directory_name=test_doc_path, extension=".txt", key="TXT")
        self.assertIsNotNone(project._input_data[0].get("TXT"))
        self.assertEqual(len(project._input_data), 1)
        # testing json
        project._input_data = []
        project.bulk_input_directory(directory_name=test_doc_path, extension=".json", key="JSON")
        self.assertIsNotNone(project._input_data[0].get("JSON"))
        self.assertEqual(len(project._input_data), 1)

    def test_bulk_input_directory_not_supported_input(self):
        project = mlp.MachineLearningProject()
        with self.assertRaises(ValueError) as context:
            project.bulk_input_directory(directory_name=files("tests") / "test_documents", extension=".doc", key=None)
            self.assertIn(".doc - extension is not supported.", context.exception)
