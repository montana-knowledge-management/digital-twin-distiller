import unittest

import numpy as np

from digital_twin_distiller.document_vectorizer import DocumentVectorizer
from importlib_resources import files

exmple_text = ["teszt szöveg első fele.", "második teszt dokumentum ami az első fele folytatása."]

example_test = ["teszt elem amit még nem látott a modell"]


class DummyFasttextModel:
    def __init__(self):
        self.word_vectors = {'teszt': np.array([-0.31889972, -0.32168077, -0.43845435, -0.05013237, 0.6502505]),
                             'szöveg': np.array([-0.24492963, -0.17807408, -0.23365464, -0.21707366, -0.33084135]),
                             'első': np.array([-0.0821762, -0.1395661, 0.29812188, 0.71711501, 0.00625835]),
                             'fele': np.array([-0.40292474, 0.10488187, 0.82989296, -0.39019406, 0.17262902]),
                             'második': np.array([-0.21978395, -0.19596792, 0.08622615, 0.19335617, -0.12503032]),
                             'dokumentum': np.array([-0.18175753, -0.22389043, -0.20449593, -0.12305195, -0.2462691]),
                             'ami': np.array([-0.17623963, 1.1600889, -0.2423223, 0.06445429, 0.0164192]),
                             'az': np.array([1.84987635, -0.01960656, 0.03841505, -0.086517, 0.04585048]),
                             'folytatása': np.array([-0.22316495, -0.18618491, -0.13372884, -0.10795644, -0.18926679])}

    def get_word_vector(self, word):
        return self.word_vectors.get(word)


class DocumentVectorizerTestCase(unittest.TestCase):
    def test_average(self):
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model
        vectorizer.fasttext_model = DummyFasttextModel()
        # creating vocabulary
        vectorizer.build_vocab(exmple_text)
        # vectorizer.build_vectors_dict(mode="fasttext")
        # creating tokenized text
        example = [txt.replace(".", "") for txt in exmple_text[0].split() if txt]
        # getting document vector
        document_vector = vectorizer.run(example, mode="average")
        print(document_vector)

        self.assertAlmostEqual(-0.2622325725, document_vector[0])
        self.assertAlmostEqual(-0.13360977, document_vector[1])

    def test_idf_weighted(self):
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model for testing
        vectorizer.fasttext_model = DummyFasttextModel()
        vectorizer.build_vocab(exmple_text)
        print(vectorizer.vocabulary)
        # vectorizer.build_vectors_dict(mode="fasttext")
        example = [txt.replace(".", "") for txt in exmple_text[0].split() if txt]
        document_vector = vectorizer.run(example, mode="idf_weighted")
        print(document_vector)
        self.assertNotEqual(-0.2622325725, document_vector[0])
        self.assertNotEqual(-0.13360977, document_vector[1])

    def test_cosine_similarity(self):
        vectorizer = DocumentVectorizer()
        vec_1 = np.array([0.1, -0.1])
        vec_2 = np.array([0.1, -0.1])
        self.assertAlmostEqual(vectorizer.cosine_similarity(vec_1, vec_2), 1.0)
        self.assertAlmostEqual(vectorizer.cosine_similarity(vec_1, vec_2 * -1), -1.0)
        vec_1 = np.array([0.1, -0.0])
        vec_2 = np.array([0.0, -0.1])
        self.assertAlmostEqual(vectorizer.cosine_similarity(vec_1, vec_2), 0.0)

    def test_docvec(self):
        vectorizer = DocumentVectorizer()
        trained_model = vectorizer.train_doc2vec_model(corpus=exmple_text, vector_dim=50, epochs=100, min_count=1)
        doc2vecvector = trained_model.infer_vector(["teszt elem amit még nem látott a modell"])
        self.assertEqual(50, len(doc2vecvector))


if __name__ == '__main__':
    unittest.main()
