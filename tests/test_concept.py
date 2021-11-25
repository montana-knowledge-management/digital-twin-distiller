import socket
import unittest

import digital_twin_distiller.concept as cpt
import digital_twin_distiller.text_readers as r
from importlib_resources import files


class TestConcepts(unittest.TestCase):
    # def test_dummy_project_ip_helper(self):
    #     # tests the basic functions of an empty abstractproject
    #     Project = cpt.AbstractProject()
    #
    #     # init
    #     self.assertIsInstance(Project, cpt.AbstractProject)
    #
    #     # check ip address helper
    #     s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    #     s.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
    #     s.connect(("<broadcast>", 0))
    #     ip = s.getsockname()[0]
    #
    #     self.assertEqual(Project.get_ip()[0], ip)

    def test_dummy_project_text_file_reader(self):
        # tests the basic functions of an empty abstractproject
        Project = cpt.AbstractProject()

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
