import unittest

import digital_twin_distiller.text_readers as rdr
import digital_twin_distiller.text_writers as wrtr


class TestWriters(unittest.TestCase):
    def test_write_json(self):
        fname = "/tmp/.sandbox.json"
        exp = {"gyümölcs": ["alma", "körte", "szilva"], "flag": True, "subkey": {"k2": 2, "k3": 3}}
        self.assertTrue(wrtr.JsonWriter().write(exp, fname))
        act = rdr.JsonReader().read(fname)
        self.assertEqual(exp, act)

    def test_write_json_bom(self):
        fname = "/tmp/.sandbox_bom.json"
        exp = {"gyümölcs": ["alma", "körte", "szilva"], "flag": True, "subkey": {"k2": 2, "k3": 3}}
        self.assertTrue(wrtr.JsonWriter().write(exp, fname, encoding="utf-8-sig"))
        act = rdr.JsonReader().read(fname)
        self.assertEqual(exp, act)

    def test_write_text(self):
        fname = "/tmp/test.txt"
        exp = "Test Text"
        self.assertTrue(wrtr.TextWriter().write(exp, fname))
        act = rdr.TextReader().read(fname)
        self.assertEqual(exp, act)
