# import unittest
# from examples.pu_learning.positive_unlabeled_learning_example import estimate_positive_count, estimate_c, \
#     compare_performances, get_potentially_positive_samples, inspect_returned_probabilities
# from sklearn.feature_extraction.text import TfidfVectorizer
# from distiller.pu_learning import PuLearning
# from sklearn.metrics import accuracy_score
#
# pu_learner = PuLearning()
# # loading exmaple dataset presented by Elkan and Noto
# pu_learner.load_example_dataset()
# # tfidf vectorization
# vectorizer = TfidfVectorizer()
# vectorizer.fit(pu_learner.P + pu_learner.U)
# p_vectors = vectorizer.transform(pu_learner.P)
# u_vectors = vectorizer.transform(pu_learner.U)
# q_vectors = vectorizer.transform(pu_learner.Q)
# pu_learner.lr_model_settings = {}
#
#
# class MyTestCase(unittest.TestCase):
#     def test_estimate_c(self):
#         pu_learner.lr_model_settings = {}
#         c_est = estimate_c(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors, model="lr")
#         self.assertAlmostEqual(c_est, 0.7372883118720045, delta=0.05)
#
#     def test_estimate_positive_count(self):
#         positive_estimate = estimate_positive_count(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors,
#                                                     random_state=1)
#         self.assertAlmostEqual(positive_estimate, 3327, delta=100)
#
#     def test_get_potentially_positive_samples(self):
#         tp, fp, fn = get_potentially_positive_samples(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors,
#                                                       model="lr")
#         self.assertAlmostEqual(tp, 269, delta=5)
#         self.assertAlmostEqual(fp, 76, delta=5)
#         self.assertAlmostEqual(fn, 79, delta=5)
#
#     def test_compare_performances(self):
#         y, yhat_orig, yhat_pu, yhat_best = compare_performances(vectorizer=vectorizer, pu_learner=pu_learner,
#                                                                 p_vectors=p_vectors, u_vectors=u_vectors, c=0.7373,
#                                                                 model="lr")
#         self.assertGreater(accuracy_score(y, yhat_orig), 0.94)
#         self.assertGreater(accuracy_score(y, yhat_pu), 0.96)
#         self.assertGreater(accuracy_score(y, yhat_best), 0.965)
#
#     def test_inspect_probabilities(self):
#         pos_prob, unlab_prob = inspect_returned_probabilities(pu_learner=pu_learner, p_vectors=p_vectors,
#                                                               u_vectors=u_vectors, model="lr")
#         self.assertEqual(len(pos_prob), 2453)
#         self.assertEqual(len(unlab_prob), 4906)
#
#
# if __name__ == '__main__':
#     unittest.main()
