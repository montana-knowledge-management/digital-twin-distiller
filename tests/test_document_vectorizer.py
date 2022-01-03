import unittest

import numpy as np
from os.path import exists
from digital_twin_distiller.document_vectorizer import DocumentVectorizer
from importlib_resources import files
from os import remove
import pytest

example_text = ["teszt szöveg első fele.", "második teszt dokumentum ami az első fele folytatása."]

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


class DummyW2VModel:
    def __init__(self):
        self.wv = {'teszt': np.array([-0.31889972, -0.32168077, -0.43845435, -0.05013237, 0.6502505]),
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
    @pytest.fixture(autouse=True)
    def _pass_fixtures(self, capsys):
        self.capsys = capsys

    def test_average(self):
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model
        vectorizer.fasttext_model = DummyFasttextModel()
        # creating vocabulary
        vectorizer.build_vocab(example_text)
        # vectorizer.build_vectors_dict(mode="fasttext")
        # creating tokenized text
        example = [txt.replace(".", "") for txt in example_text[0].split() if txt]
        # getting document vector
        document_vector = vectorizer.run(example, mode="average")
        print(document_vector)

        self.assertAlmostEqual(-0.2622325725, document_vector[0])
        self.assertAlmostEqual(-0.13360977, document_vector[1])

    def test_idf_weighted(self):
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model for testing
        vectorizer.fasttext_model = DummyFasttextModel()
        vectorizer.build_vocab(example_text)
        print(vectorizer.vocabulary)
        # vectorizer.build_vectors_dict(mode="fasttext")
        example = [txt.replace(".", "") for txt in example_text[0].split() if txt]
        document_vector = vectorizer.run(example, mode="idf_weighted")
        print(document_vector)
        self.assertNotEqual(-0.2622325725, document_vector[0])
        self.assertNotEqual(-0.13360977, document_vector[1])

    def test_build_vectors_dict(self):
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model for testing
        vectorizer.fasttext_model = DummyFasttextModel()
        vectorizer.build_vocab(example_text)
        vectorizer.build_vectors_dict(mode="fasttext")
        self.assertDictEqual(vectorizer.vectors_dict, vectorizer.fasttext_model.word_vectors)

    def test_build_vectors_dict_raise(self):
        with self.assertRaises(ValueError) as context:
            vectorizer = DocumentVectorizer()
            vectorizer.build_vectors_dict(mode="nonexisting")
            self.assertEqual(context, "Wrong mode given! Please choose from 'gensim' or 'fasttext'!")

        with self.assertRaises(ValueError) as context:
            vectorizer = DocumentVectorizer()
            vectorizer.build_vectors_dict()
            self.assertEqual(context, "Missing loaded fasttext model. Please load fasttext model!")

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
        model_path_to_save = "/tmp/test.bin"
        vectorizer = DocumentVectorizer()
        trained_model = vectorizer.train_doc2vec_model(corpus=example_text, vector_dim=50, epochs=100, min_count=1,
                                                       model_path_to_save=model_path_to_save)
        self.assertTrue(exists(model_path_to_save))

        vectorizer = DocumentVectorizer()
        vectorizer.load_doc2vec_model(model_path_to_save)

        doc2vec_vector = trained_model.infer_vector(example_test)
        loaded_doc2vec_vector = vectorizer.doc2vec_model.infer_vector(example_test)
        # removing trained model
        remove(model_path_to_save)
        self.assertEqual(50, len(doc2vec_vector))
        self.assertEqual(50, len(loaded_doc2vec_vector))
        # checking whether loaded and original models are the same
        self.assertListEqual(doc2vec_vector.tolist(), loaded_doc2vec_vector.tolist())
        # getting vector form via vectorizer object
        document_vector = vectorizer.run(example_test, mode="doc2vec")
        # checking whether run method gives the same result for the same input
        self.assertListEqual(doc2vec_vector.tolist(), document_vector.tolist())

    def test_keep_n_most_similar_to_average(self):
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model
        vectorizer.fasttext_model = DummyFasttextModel()
        # creating vocabulary
        vectorizer.build_vocab(example_text)

        with self.assertRaises(ValueError) as context:
            vectorizer.keep_n_most_similar_to_average([""], mode="notimplemented")
            self.assertEqual(context,
                             'The mode is currently not implemented, please choose from "average" and "idf_weighted"!')

        example_tokenized = [txt.replace(".", "") for txt in example_text[0].split() if txt]
        mean_of_most_similar = vectorizer.keep_n_most_similar_to_average(example_tokenized, nr_of_words_to_keep=2,
                                                                         mode="average")
        szoveg_elso_vec_lst = np.array(
            [vectorizer.fasttext_model.word_vectors.get("szöveg"), vectorizer.fasttext_model.word_vectors.get("első")])
        # asserting that the words szöveg and első are the closest two to the mean vector
        self.assertListEqual(mean_of_most_similar.tolist(), szoveg_elso_vec_lst.mean(axis=0).tolist())

    def test_keep_n_most_similar_to_idf(self):
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model
        vectorizer.fasttext_model = DummyFasttextModel()
        # creating vocabulary
        vectorizer.build_vocab(example_text)

        example_tokenized = [txt.replace(".", "") for txt in example_text[0].split() if txt]
        mean_of_most_similar = vectorizer.keep_n_most_similar_to_average(example_tokenized, nr_of_words_to_keep=2,
                                                                         mode="idf_weighted")
        szoveg_elso_vec_lst = np.array(
            [vectorizer.fasttext_model.word_vectors.get("szöveg") * vectorizer.idf[vectorizer.vocabulary.get("szöveg")],
             vectorizer.fasttext_model.word_vectors.get("első") * vectorizer.idf[vectorizer.vocabulary.get("első")]])
        # asserting that the words szöveg and első are the closest two to the mean vector
        self.assertListEqual(mean_of_most_similar.tolist(), szoveg_elso_vec_lst.mean(axis=0).tolist())

        mean_of_most_similar_avg = vectorizer.keep_n_most_similar_to_average(example_tokenized, nr_of_words_to_keep=2,
                                                                             mode="average")
        # checking if idf vectors and avg vectors are different
        self.assertNotEqual(mean_of_most_similar_avg.mean(), mean_of_most_similar.mean())

    def test_save_and_load_vectors_dict(self):
        path_to_save = "/tmp/vectors.json"
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model
        vectorizer.fasttext_model = DummyFasttextModel()
        # creating vocabulary
        vectorizer.build_vocab(example_text)
        vectorizer.build_vectors_dict(mode="fasttext")
        # saving vectors dict
        vectorizer.save_vectors_dictionary(path_to_save)
        self.assertTrue(exists(path_to_save))
        print(vectorizer.vectors_dict)
        # loading vectors dict
        vectorizer.vectors_dict = None
        vectorizer.load_vectors_dictionary(path_to_save)
        print(vectorizer.vectors_dict)
        self.assertEqual(vectorizer.vectors_dict.get("szöveg").tolist(),
                         vectorizer.fasttext_model.word_vectors.get("szöveg").tolist())
        # cleaning up
        remove(path_to_save)
        self.assertFalse(exists(path_to_save))

    def test_calculate_doc2vec(self):
        with self.assertRaises(ValueError) as context:
            vectorizer = DocumentVectorizer()
            vectorizer.calculate_doc2vec([""])
            self.assertEqual(context, "Missing Doc2Vec model!")

    def test_gensim_missing(self):
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model
        vectorizer.fasttext_model = DummyFasttextModel()
        # creating vocabulary
        vectorizer.build_vocab(example_text)
        with self.assertRaises(ValueError) as context:
            vectorizer.build_vectors_dict(mode="gensim")
            self.assertEqual(context, "Missing loaded gensim model. Please load gensim model!")

    def test_gensim_model(self):
        vectorizer = DocumentVectorizer()
        # loading dummy fasttext model
        vectorizer.gensim_model = DummyW2VModel()
        # creating vocabulary
        vectorizer.build_vocab(example_text)
        vectorizer.build_vectors_dict(mode="gensim")
        self.assertDictEqual(vectorizer.vectors_dict, vectorizer.gensim_model.wv)

        vocab = list(vectorizer.vocabulary)
        vocab.append("pentakosziomedimnosz")
        vectorizer.vocabulary = set(vocab)
        vectorizer.build_vectors_dict(mode="gensim")
        captured = self.capsys.readouterr()
        self.assertEqual(captured.out.replace("\n", ""),
                         "The word pentakosziomedimnosz was missing from the gensim model, not putting into vectors dict.")


if __name__ == '__main__':
    unittest.main()
