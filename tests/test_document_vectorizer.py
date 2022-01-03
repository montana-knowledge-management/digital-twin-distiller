import unittest

import numpy as np

from digital_twin_distiller.document_vectorizer import DocumentVectorizer
from importlib_resources import files

fasttext_path = "/media/csanyig/C8CC19CCCC19B622/Users/csanyig/Montana/NLP/Word_embeddings/cc.hu.2.bin"

exmple_text = ["teszt szöveg első fele.", "második teszt dokumentum ami az első fele folytatása."]

example_test = ["teszt elem amit még nem látott a modell"]


class MyTestCase(unittest.TestCase):
    def test_average(self):
        vectorizer = DocumentVectorizer()
        vectorizer.load_fasttext_model(fasttext_path)
        # vectorizer.build_vocab(exmple_text)
        # vectorizer.build_vectors_dict(mode="fasttext")
        example = [txt for txt in exmple_text[0].split() if txt]
        document_vector = vectorizer.run(example, mode="average")
        print(document_vector)
        self.assertAlmostEqual(-0.04909943, document_vector[0])
        self.assertAlmostEqual(-0.07539089, document_vector[1])

    def test_idf_weighted(self):
        vectorizer = DocumentVectorizer()
        vectorizer.load_fasttext_model(fasttext_path)
        vectorizer.build_vocab(exmple_text)
        # vectorizer.build_vectors_dict(mode="fasttext")
        example = [txt for txt in exmple_text[0].split() if txt]
        document_vector = vectorizer.run(example, mode="idf_weighted")
        print(document_vector)
        self.assertNotEqual(-0.04909943, document_vector[0])
        self.assertNotEqual(-0.07539089, document_vector[1])

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
