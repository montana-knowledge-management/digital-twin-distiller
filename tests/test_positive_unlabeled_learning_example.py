import unittest
from examples.text_processing.pu_learning.positive_unlabeled_learning_example import estimate_positive_count, \
    estimate_c, \
    compare_performances, get_potentially_positive_samples, inspect_returned_probabilities
from sklearn.feature_extraction.text import TfidfVectorizer
from digital_twin_distiller.pu_learning import PuLearning
from sklearn.metrics import accuracy_score

pu_learner = PuLearning()
pu_learner.real_c = 0.6
# loading exmaple dataset presented by twenty newsgroups dataset
pu_learner.load_example_dataset()
# tfidf vectorization
vectorizer = TfidfVectorizer()
vectorizer.fit(pu_learner.P + pu_learner.U)
p_vectors = vectorizer.transform(pu_learner.P)
u_vectors = vectorizer.transform(pu_learner.U)
q_vectors = vectorizer.transform(pu_learner.Q)
pu_learner.lr_model_settings = {"class_weight": "balanced"}


class MyTestCase(unittest.TestCase):
    def test_estimate_c(self):
        c_est = estimate_c(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors, model="lr")
        self.assertAlmostEqual(c_est, pu_learner.real_c, delta=0.05)  # Elkan Noto: 0.7372883118720045

    def test_estimate_positive_count(self):
        positive_estimate = estimate_positive_count(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors,
                                                    random_state=1)
        self.assertAlmostEqual(positive_estimate, len(pu_learner.P) + len(pu_learner.Q), delta=100)  # Elkan Noto: 3327

    def test_get_potentially_positive_samples(self):
        tp, fp, fn = get_potentially_positive_samples(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors,
                                                      model="lr")
        self.assertAlmostEqual(tp, 343, delta=5)  # 269
        self.assertAlmostEqual(fp, 128, delta=5)  # 76
        self.assertAlmostEqual(fn, 53, delta=5)  # 79

    def test_compare_performances(self):
        y, yhat_orig, yhat_pu, yhat_best = compare_performances(vectorizer=vectorizer, pu_learner=pu_learner,
                                                                p_vectors=p_vectors, u_vectors=u_vectors, c=0.6,
                                                                # c=0.7373,
                                                                model="lr")
        self.assertGreater(accuracy_score(y, yhat_orig), 0.84)  # 0.94
        self.assertGreater(accuracy_score(y, yhat_pu), 0.90)  # 0.96
        self.assertGreater(accuracy_score(y, yhat_best), 0.99)  # 0.965

    def test_inspect_probabilities(self):
        pos_prob, unlab_prob = inspect_returned_probabilities(pu_learner=pu_learner, p_vectors=p_vectors,
                                                              u_vectors=u_vectors, model="lr")
        self.assertEqual(len(pos_prob), 594)  # 2453
        self.assertEqual(len(unlab_prob), 1392)  # 4906


if __name__ == '__main__':
    unittest.main()
