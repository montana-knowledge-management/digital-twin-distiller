import unittest

import numpy as np
import scipy.sparse.csr

from digital_twin_distiller.pu_learning import PuLearning
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import MultinomialNB
from sklearn.svm import SVC, LinearSVC


class PuLearningTestCase(unittest.TestCase):
    # def test_c_estimation(self):
    #     difference = 0.15
    #     pu_learner = PuLearning()
    #     # loading P U dataset from Elkan and Noto
    #     pu_learner.load_example_dataset()
    #     vectorizer = TfidfVectorizer()
    #     vectorizer.fit(pu_learner.P + pu_learner.U)
    #     p_vectors = vectorizer.transform(pu_learner.P)
    #     u_vectors = vectorizer.transform(pu_learner.U)
    #     pu_learner.estimate_c(P=p_vectors, U=u_vectors, model="lr", n_splits=10, random_state=1)
    #     print(pu_learner.real_c, pu_learner.c)
    #     self.assertLess(pu_learner.real_c - pu_learner.c, difference)

    def test_estimate_positive_count(self):
        c = 0.75
        count = 1000
        P = [1] * count
        pu_learner = PuLearning()
        pu_learner.c = c
        self.assertEqual(pu_learner.estimate_positive_count(P, []), 1333)

    def test_choose_model(self):
        random_state = 7
        pu_learner = PuLearning()
        model1 = LogisticRegression()
        model2 = "lr"
        model3 = "svc"
        model4 = MultinomialNB()
        self.assertTrue(isinstance(pu_learner._choose_model(model=model1), LogisticRegression))
        self.assertTrue(isinstance(pu_learner._choose_model(model=model2), LogisticRegression))
        self.assertTrue(isinstance(pu_learner._choose_model(model=model3), SVC))
        self.assertTrue(isinstance(pu_learner._choose_model(model=model4), MultinomialNB))
        self.assertEqual(pu_learner._choose_model(random_state=random_state).random_state, random_state)

    def test_choose_model_raising_error(self):
        pu_learner = PuLearning()
        # wrong model shorthand
        with self.assertRaises(ValueError) as context:
            pu_learner._choose_model(model="not stated")
        self.assertTrue("Wrong model type given. Please choose from ['lr','svc'] or use a model object where the 'fit' "
                        "and 'predict_proba' functions are implemented!" in str(context.exception))
        # only fit method is implemented in case of Linear SVC
        with self.assertRaises(NotImplementedError) as context:
            pu_learner._choose_model(model=LinearSVC())
        self.assertTrue(
            "Please use a model where both 'fit' and 'predict_proba' functions are implemented!" in str(
                context.exception))

    # def test_load_example_dataset(self):
    #     pu_learner = PuLearning()
    #     pu_learner.load_example_dataset()
    #     self.assertEqual(len(pu_learner.P), 2453)
    #     self.assertEqual(len(pu_learner.Q), 348)
    #     self.assertEqual(len(pu_learner.N), 4558)
    #     self.assertEqual(len(pu_learner.U), 4906)
    #     self.assertAlmostEqual(pu_learner.real_c, 0.8758, delta=10e-04)
    #     self.assertEqual(pu_learner.c, 0.0)

    def test_create_training_data_list(self):
        pu_learner = PuLearning()
        P = [1, 2, 3]
        U = [5, 5, 5, 2]
        X, y = pu_learner.create_training_data(P, U)
        self.assertListEqual(X, [1, 2, 3, 5, 5, 5, 2])
        self.assertListEqual(y.tolist(), [1, 1, 1, 0, 0, 0, 0])
        self.assertTrue(isinstance(X, list))
        self.assertTrue(isinstance(y, np.ndarray))

    def test_create_training_data_exception(self):
        pu_learner = PuLearning()
        P = (1, 2, 3)
        U = [5, 5, 5, 2]
        with self.assertRaises(ValueError) as context:
            pu_learner.create_training_data(P, U)
        self.assertTrue('P and U are not the following types: list, numpy.ndarray' in str(context.exception))

    # def test_create_training_data_sparse(self):
    #     pu_learner = PuLearning()
    #     pu_learner.load_example_dataset()
    #     vectorizer = TfidfVectorizer()
    #     vectorizer.fit(pu_learner.P + pu_learner.U)
    #     p_vectors = vectorizer.transform(pu_learner.P)
    #     u_vectors = vectorizer.transform(pu_learner.U)
    #     X, y = pu_learner.create_training_data(p_vectors, u_vectors)
    #     self.assertTrue(isinstance(p_vectors, scipy.sparse.csr.csr_matrix))
    #     self.assertEqual(X.shape, (7359, 44575))
    #     self.assertTrue(isinstance(y, np.ndarray))
    #     self.assertEqual(np.sum(y), p_vectors.shape[0])

    def test_create_training_data_numpy(self):
        pu_learner = PuLearning()
        P = np.array([1, 2, 3])
        U = np.array([5, 5, 5, 2])
        X, y = pu_learner.create_training_data(P, U)
        self.assertTrue(np.all(X == np.array([1, 2, 3, 5, 5, 5, 2])))
        self.assertListEqual(y.tolist(), [1, 1, 1, 0, 0, 0, 0])

    def test_init(self):
        pu_learner = PuLearning()
        self.assertEqual(pu_learner.c, 0.0)
        self.assertEqual(pu_learner.P, [])
        self.assertEqual(pu_learner.Q, [])
        self.assertEqual(pu_learner.N, [])
        self.assertEqual(pu_learner.U, [])
        self.assertDictEqual(pu_learner.svc_model_settings, {"kernel": "linear", "probability": True, "C": 1})
        self.assertDictEqual(pu_learner.lr_model_settings, {'class_weight': 'balanced'})

    # def test_train_test_split(self):
    #     pu_learner = PuLearning()
    #     pu_learner.load_example_dataset()
    #     vectorizer = TfidfVectorizer()
    #     vectorizer.fit(pu_learner.P + pu_learner.U)
    #     p_vectors = vectorizer.transform(pu_learner.P)
    #     u_vectors = vectorizer.transform(pu_learner.U)
    #     X, y = pu_learner.create_training_data(p_vectors, u_vectors)
    #     X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True, stratify=y, random_state=1,
    #                                                         test_size=0.2)
    #     model = LogisticRegression(random_state=1)
    #     c_estimate = pu_learner._train_classifier_to_estimate_c(X_train, y_train, X_test, y_test, model)
    #     self.assertAlmostEqual(c_estimate, 0.7406, delta=5e-02)

    def test_train_estimator(self):
        pu_learner = PuLearning()
        P = np.array([[1], [2], [3]])
        U = np.array([[5], [5], [5], [2]])
        model = pu_learner.train_estimator(P, U, model="lr", random_state=1)
        result = model.predict_proba(np.array([[1.5]]))
        self.assertEqual(model.random_state, 1)
        self.assertAlmostEqual(result[0][1], 0.798656889234633, delta=0.01)

        lr_model = LogisticRegression(random_state=1, class_weight="balanced")
        X, y = pu_learner.create_training_data(P, U)
        lr_model.fit(X, y)
        self.assertAlmostEqual(result[0][1], lr_model.predict_proba(np.array([[1.5]]))[0][1], delta=0.01)

    def test_get_potential_positives(self):
        pu_learner = PuLearning()
        P = np.array([[1], [2], [3]])
        U = np.array([[5], [5], [5], [2]])
        potential_positive_index = pu_learner.get_potential_positives(P=P, U=U, model="lr", n_splits=2,
                                                                      random_state=1)
        self.assertEqual(potential_positive_index, [3])
        self.assertEqual(U[potential_positive_index], [2])
        self.assertEqual(pu_learner.estimate_positive_count(P=P, U=U, n_splits=2, random_state=1), 4)


if __name__ == '__main__':
    unittest.main()
