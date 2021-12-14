import unittest

from importlib_resources import files

from digital_twin_distiller import text_readers


class TestReaders(unittest.TestCase):
    def test_read_pdf(self):
        reader = text_readers.PdfReader()
        act = reader.read(file_name=str(files("tests") / "test_documents" / "Moore.pdf"))
        self.assertEqual(reader.ok, act["status"])
        self.assertIn(reader.text_key, set(act.keys()))
        self.assertIn(reader.status, set(act.keys()))
        self.assertEqual("Law of Moore", act.get("title"))
        self.assertIn("While Moore did not use", act["Text"])
        self.assertIn(
            "microprocessor prices, the increase in memory capacity (RAM and flash), the improvement of",
            act["Text"].replace("\n", ""),
        )
