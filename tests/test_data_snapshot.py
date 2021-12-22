import unittest
from os import path
from os import remove
from zipfile import ZipFile

from digital_twin_distiller.data_store.data_snapshot import DataSnapshot


class TestDataSnapshot(unittest.TestCase):
    def test_basic_save_function(self):
        data = [{"dict_1": "Hello!", "FileName": "Name"}, {"dict_2": "Hello!!"}]
        DataSnapshot.save(data, "temp")
        # is the file is created
        self.assertTrue(path.isfile("temp.zip"))

        with ZipFile("temp.zip", "r") as zipObj:
            listOfiles = zipObj.namelist()
            self.assertIn("Name.json", listOfiles)
            self.assertEqual(2, len(listOfiles))
        try:
            remove("temp.zip")
        except OSError:
            print("The test file didn't created successfully!")

    def test_basic_load_function(self):
        data = [{"dict_1": "Hello!", "FileName": "Name"}, {"dict_2": "Hello!!"}]
        DataSnapshot.save(data, "temp")
        # is the file is created

        stack = DataSnapshot.load_stack("temp.zip")
        print(stack)

        try:
            remove("temp.zip")
        except OSError:
            print("The test file didn't created successfully!")
