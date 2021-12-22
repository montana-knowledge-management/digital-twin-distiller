"""
Positive and unlabeled learning based on Elkan and Noto: Learning Classifiers from Only Positive and Unlabeled Data
(https://cseweb.ucsd.edu/~elkan/posonly.pdf) paper.
The dataset is available from http://cseweb.ucsd.edu/~elkan/posonly/.
However, the examples were modified to use the twenty newsgroups dataset.
"""
import numpy as np
from digital_twin_distiller.pu_learning import PuLearning
from sklearn.feature_extraction.text import TfidfVectorizer
import scipy
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score


def compare_performances(vectorizer, pu_learner, p_vectors, u_vectors, c, model="lr"):
    """Example code for the effect of estimation of c for the P vs U, P vs U wighted and ideal P+Q vs N"""
    pu_learner.c = c
    pu_model = pu_learner.train_estimator(P=p_vectors, U=u_vectors, model=model, random_state=1)
    threshold = 0.5 * pu_learner.c
    real_positive_vectors = vectorizer.transform(pu_learner.P + pu_learner.Q)
    real_negative_vectors = vectorizer.transform(pu_learner.N)
    non_pu_model = pu_learner.train_estimator(P=real_positive_vectors, U=real_negative_vectors, model=model,
                                              random_state=1)
    X = scipy.sparse.vstack((real_positive_vectors, real_negative_vectors))
    y = [1] * real_positive_vectors.shape[0] + [0] * real_negative_vectors.shape[0]
    yhat_pu = [prob[1] >= threshold for prob in pu_model.predict_proba(X)]
    yhat_orig = [prob[1] >= 0.5 for prob in pu_model.predict_proba(X)]
    yhat_best = [prob[1] >= 0.5 for prob in non_pu_model.predict_proba(X)]
    print(
        "\n======================\nResults after training on Positive and unlabeled as negatives (P vs. Q+N), with threshold 0.5"
        "\n======================\n")
    print(classification_report(y_true=y, y_pred=yhat_orig))
    print("\nConfusion matrix:\n")
    print(confusion_matrix(y_true=y, y_pred=yhat_orig))
    print("\nAccuracy: {:.2f}%\n".format(accuracy_score(y, yhat_orig) * 100))
    print(
        "\n======================\nResults after training on Positive and unlabeled as negatives (P vs. Q+N), with threshold {}"
        "\n======================\n".format(threshold))
    print(classification_report(y_true=y, y_pred=yhat_pu))
    print("\nConfusion matrix:\n")
    print(confusion_matrix(y_true=y, y_pred=yhat_pu))
    print("\nAccuracy: {:.2f}%\n".format(accuracy_score(y, yhat_pu) * 100))
    print(
        "\n======================\nResults after training on all positive and all negatives (P+Q vs. N)\n======================\n")
    print(classification_report(y_true=y, y_pred=yhat_best))
    print("\nConfusion matrix:\n")
    print(confusion_matrix(y_true=y, y_pred=yhat_best))
    print("\nAccuracy: {:.2f}%\n".format(accuracy_score(y, yhat_best) * 100))
    return y, yhat_orig, yhat_pu, yhat_best


def estimate_c(pu_learner, p_vectors, u_vectors, model="lr"):
    """Example code how to estimate value c"""
    pu_learner.estimate_c(P=p_vectors, U=u_vectors, model=model, n_splits=10, random_state=1)
    print("\n=================\nReal value of c: {}, estimated value of c: {}\n=================\n".format(
        pu_learner.real_c, pu_learner.c))
    return pu_learner.c


def estimate_positive_count(pu_learner, p_vectors, u_vectors, random_state=None):
    """Example code for estimating the number of positive samples in U dataset"""
    estimated_positive_count = pu_learner.estimate_positive_count(P=p_vectors, U=u_vectors, random_state=random_state)
    print(
        "\n=================\n"
        "Estimated number of positive samples: {}, real number of positive samples: {}, difference: {}"
        "\n=================\n".format(estimated_positive_count, len(pu_learner.P) + len(pu_learner.Q),
                                       estimated_positive_count - (len(pu_learner.P) + len(pu_learner.Q))))
    return estimated_positive_count


def get_potentially_positive_samples(pu_learner, p_vectors, u_vectors, model="lr"):
    """Example code for extracting the potential positives from the unlabeled dataset"""
    potential_indices = pu_learner.get_potential_positives(P=p_vectors, U=u_vectors, model=model,
                                                           random_state=1, plot=False)
    good_finds = [idx < len(pu_learner.Q) for idx in potential_indices]
    tp = sum(good_finds)
    fp = len(good_finds) - sum(good_finds)
    fn = len(pu_learner.Q) - sum(good_finds)
    print("\n=================\n"
          "How well the method could find the positive samples from the Unlabeled dataset:"
          "\n=================\n")
    print("Correctly found elements (TP): {}, wrong positives (FP): {}, missed examples (FN): {}\n"
          "Precision: {:.2f}%, Recall: {:.2f}%".format(tp, fp, fn, 100 * tp / (tp + fp), 100 * tp / (tp + fn)))
    return tp, fp, fn


def inspect_returned_probabilities(pu_learner, p_vectors, u_vectors, model="lr"):
    """Example code for utilizing the returned probabilities after training a non-traditional classifier on P and U"""
    p_probs, u_probs = pu_learner.train_and_reveal_probabilities(P=p_vectors, U=u_vectors, model=model,
                                                                 random_state=1)
    p_probs = np.array(p_probs)
    u_probs = np.array(u_probs)
    print("\n=================\nPositives: Avg: {:.2f}%, Std: {:.2f}%, Min: {:.2f}%, Max: {:.2f}%".format(
        np.mean(p_probs) * 100,
        np.std(p_probs) * 100,
        np.min(p_probs) * 100,
        np.max(p_probs) * 100))
    print("Unlabeled: Avg: {:.2f}%, Std: {:.2f}%, Min: {:.2f}%, Max: {:.2f}%".format(np.mean(u_probs) * 100,
                                                                                     np.std(u_probs) * 100,
                                                                                     np.min(u_probs) * 100,
                                                                                     np.max(u_probs) * 100))
    return p_probs, u_probs


if __name__ == "__main__":
    pu_learner = PuLearning()
    # choose from "lr" or "svc" or use model where fit and predic_proba methods are available
    model = "lr"
    # loading exmaple dataset presented by Elkan and Noto
    pu_learner.load_example_dataset()
    pu_learner.lr_model_settings.update({"class_weight": "balanced", "n_jobs": -1})
    pu_learner.svc_model_settings.update({"class_weight": "balanced"})
    # vectorizing
    vectorizer = TfidfVectorizer()
    vectorizer.fit(pu_learner.P + pu_learner.U)
    p_vectors = vectorizer.transform(pu_learner.P)
    u_vectors = vectorizer.transform(pu_learner.U)
    q_vectors = vectorizer.transform(pu_learner.Q)

    estimate_c(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors, model=model)
    estimate_positive_count(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors, random_state=1)

    # compare performances of the given model ideal training, training on PU sets with 0.5 threshold and with 0.5c threshold
    compare_performances(vectorizer, pu_learner, p_vectors, u_vectors, pu_learner.c, model=model)

    get_potentially_positive_samples(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors,
                                     model=model)

    inspect_returned_probabilities(pu_learner=pu_learner, p_vectors=p_vectors, u_vectors=u_vectors,
                                   model=model)
