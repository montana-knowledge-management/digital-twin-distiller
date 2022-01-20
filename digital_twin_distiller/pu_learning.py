"""
Positive and unlabeled learning based on Elkan and Noto: Learning Classifiers from Only Positive and Unlabeled Data
(https://cseweb.ucsd.edu/~elkan/posonly.pdf) paper.
The dataset is available from http://cseweb.ucsd.edu/~elkan/posonly/.
However, the examples were modified to use the twenty newsgroups dataset.
"""
from digital_twin_distiller.ml_project import AbstractTask
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from importlib_resources import files
from digital_twin_distiller.data_store.data_snapshot import DataSnapshot
import os
import scipy
from sklearn.model_selection import StratifiedKFold
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import fetch_20newsgroups


class PuLearning(AbstractTask):
    c = 0.0
    # known positive samples
    P = []
    # positive samples in unlabeled
    Q = []
    # negatives in unlabeled
    N = []
    # U = N union Q
    U = []
    # settings used in SVC model
    svc_model_settings = {"kernel": "linear", "probability": True, "C": 1}
    # settings used in Logistic Regression model
    lr_model_settings = {"class_weight": "balanced"}

    def _choose_model(self, model="lr", random_state=None):
        """
        Method for handling built in and custom models for analysis.
        :param model: Type of model used during the estimation:
                            lr: Logistic Regression
                            svc: linear kernel support vector classifier
                            custom model where 'fit' and 'predict_proba' methods are implemented
        :param random_state: random state used during training
        :return: model object that implements both fit and predict_proba methods
        """
        if model not in ["lr", "svc"] and not isinstance(model, str):

            if hasattr(model, 'fit') and hasattr(model, 'predict_proba'):
                if hasattr(model, 'random_state'):
                    model.random_state = random_state
                custom_model = model
            else:
                raise NotImplementedError(
                    "Please use a model where both 'fit' and 'predict_proba' functions are implemented!")
        elif model == "lr":
            custom_model = LogisticRegression(**self.lr_model_settings, random_state=random_state)
        elif model == "svc":
            custom_model = SVC(**self.svc_model_settings, random_state=random_state)
        else:
            raise ValueError("Wrong model type given. Please choose from ['lr','svc'] or use a model object where the "
                             "'fit' and 'predict_proba' functions are implemented!")
        return custom_model

    # def load_example_dataset(self):
    #     """
    #     Loads example dataset from the paper of Elkan and Noto paper cited in the headline.
    #     :return:
    #     """
    #     pu_path = files("distiller") / "resources" / "pu_learning"
    #     q_path = os.path.join(pu_path, "Q.zip")
    #     p_path = os.path.join(pu_path, "P.zip")
    #     n_path = os.path.join(pu_path, "N.zip")
    #     self.Q = DataSnapshot.load_stack(q_path)
    #     self.P = DataSnapshot.load_stack(p_path)
    #     self.N = DataSnapshot.load_stack(n_path)
    #     self.U = self.Q + self.N
    #     self.real_c = len(self.P) / (len(self.P) + len(self.Q))

    def load_example_dataset(self):
        """
        Loads example dataset from sklearn 20 newsgroups dataset.
        :return:
        """
        self.P = fetch_20newsgroups(subset='train', categories=["rec.autos"]).data
        self.Q = fetch_20newsgroups(subset='test', categories=["rec.autos"]).data
        self.N = fetch_20newsgroups(subset="all", categories=["rec.motorcycles"]).data
        self.U = self.Q + self.N
        self.real_c = len(self.P) / (len(self.P) + len(self.Q))

    def create_training_data(self, P, U):
        """
        Modifies the positive and unlabeled datasets into a machine learning format.
        :param P: positive samples, either list or sparse
        :param U: unlabeled samples, either list or sparse
        :return: X: P and U sets merged together
                 y: numpy array of labels, 1 for instances of P, 0 otherwise
        """
        if isinstance(P, list) and isinstance(U, list):
            X = P + U
            y = [1] * len(P) + [0] * len(U)
        elif isinstance(P, scipy.sparse.csr.csr_matrix) and isinstance(U, scipy.sparse.csr.csr_matrix):
            X = scipy.sparse.vstack((P, U))
            y = [1] * P.shape[0] + [0] * U.shape[0]
        elif isinstance(P, np.ndarray) and isinstance(U, np.ndarray):
            X = np.concatenate((P, U), axis=0)
            y = [1] * len(P) + [0] * len(U)
        else:
            raise ValueError("P and U are not the following types: list, numpy.ndarray, scipy.sparse.csr.csr_matrix,"
                             " or P and U are not the same types.")
        y = np.array(y)
        return X, y

    def train_estimator(self, P: list, U: list, model="lr", random_state=None):
        """
        Trains an estimator on the positive and the unlabeled dataset by counting unlabeled data as negatives.
        :param P: positive samples, either list or sparse
        :param U: unlabeled samples, either list or sparse
        :param model: Type of model used during the estimation:
                            lr: Logistic Regression
                            svc: linear kernel support vector classifier
                            custom model where 'fit' and 'predict_proba' methods are implemented
        :param random_state: random state used during training
        :return:
        """
        X, y = self.create_training_data(P=P, U=U)
        chosen_model = self._choose_model(model=model, random_state=random_state)
        chosen_model.fit(X, y)
        return chosen_model

    def train_and_reveal_probabilities(self, P: list, U: list, model="lr", random_state=None):
        """
        Trains a classifier on the P and U datasets as positives and negatives.
        :param P: positive samples in vectorized format
        :param U: unlabeled samples in vectorized format
        :param model: Type of model used during the estimation:
                            lr: Logistic Regression
                            svc: linear kernel support vector classifier
                            custom model where 'fit' and 'predict_proba' methods are implemented
        :param random_state: random state to be used during the K-Fold cross validation and ML model training
        :return: predicted probabilities for P set, predicted probabilities for U set
        """
        model = self.train_estimator(P, U, model=model, random_state=random_state)
        predicted_probs_p = [prob[1] for prob in model.predict_proba(P)]
        predicted_probs_u = [prob[1] for prob in model.predict_proba(U)]

        return predicted_probs_p, predicted_probs_u

    def get_potential_positives(self, P: list, U: list, threshold=None, model="lr", n_splits=10,
                                random_state=None, plot=False):
        """
        Returns the indices if potentially positive data from the U dataset. To achieve this the threshold 0.5c is used
        where c is the estimated probability of a positive sample being labeled. If threshold is set manually, it is
        used as it is as a threshold as the minimum probability for a positive label.
        Otherwise if c parameter has been calculated earlier and available inside the class is used. In case this is 0,
        c is esimated using the P and U datasets.
        :param P: positive samples in vectorized format
        :param U: unlabeled samples in vectorized format
        :param threshold: threshold classification [0,1], if not given and c has not been estimated, it is estimated
                          automatically based on the input data, and the classification threshold is set to 0.5c.
        :param model: Type of model used during the estimation:
                            lr: Logistic Regression
                            svc: linear kernel support vector classifier
                            custom model where 'fit' and 'predict_proba' methods are implemented
        :param random_state: random state to be used during the K-Fold cross validation and ML model training
        :param n_splits: int, default=10, number os splits used for cross validation during the estimation
        :param plot: boolean,
        :return: indices of the potential positive samples
        """
        # trainin a classifier on the P and U datasets, and getting the probabilities
        predicted_probs_p, predicted_probs_u = self.train_and_reveal_probabilities(P, U, model, random_state)
        if not threshold:
            if self.c == 0:
                self.estimate_c(P=P, U=U, n_splits=n_splits, random_state=random_state, model=model)
            # setting classification threshold to 0.5c
            threshold = self.c * 0.5
        # Plotting two histograms
        if plot:
            plt.hist(predicted_probs_p, bins=50, range=(0, 1), color="red", alpha=0.5, label="Members of P dataset")
            plt.hist(predicted_probs_u, bins=50, range=(0, 1), color="blue", alpha=0.5, label="Members of U dataset")
            plt.hist(predicted_probs_u, bins=50, range=(0, 1), color="blue", alpha=0.5, label="Members of U dataset")
            plt.axvline(x=threshold, color="green", label="Threshold")
            plt.legend(loc='upper right')
            plt.xlabel("Probabilities given by the non-traditional classifier")
            plt.ylabel("Count of documents (pcs)")
            plt.show()
        potential_good_indices = [idx for idx, prob in enumerate(predicted_probs_u) if prob > threshold]
        return potential_good_indices

    @staticmethod
    def _train_classifier_to_estimate_c(X_train, y_train, X_test, y_test, model):
        """
        Calculates the estimation of c parameter for a given machine learning model and train-teest splits.
        :param X_train: vectors of training data
        :param y_train: training labels
        :param X_test: vectors of test data
        :param y_test: test labels
        :param model: machine learning model, in the paper Logistic Regression or linear support vector machine
        :return: estimation of c
        """
        # fitting model on train data
        model.fit(X_train, y_train)
        # predicting probabilities of the test set
        predicted_probs = model.predict_proba(X_test)
        # getting the indices of known true labels
        true_indices = [idx for idx, elem in enumerate(y_test) if elem > 0]
        # collecting the probabilities of positive classification
        true_predictions = [predicted_probs[index][1] for index in true_indices]

        # avoiding division by zero
        if len(true_indices) > 0:
            c = sum(true_predictions) / len(true_indices)
            print("Averaged on {} predictions, c is: {:.2f}%.".format(len(true_indices), c * 100))
        else:
            raise ValueError("Validation set is empty!")
        return c

    def estimate_c(self, P: list, U: list, n_splits=10, random_state=None, model="lr"):
        """
        Estimates the constant probability that a positive example is labeled: c=p(s=1|y=1).
        :param P: list of vectors of positive samples
        :param U: list of vectors of unlabeled samples
        :param model_type: Type of model used during the estimation:
                            lr: Logistic Regression
                            svc: linear kernel support vector classifier
        :param n_splits: default=10, number of splits for cross validation
        :param random_state: int, random state to be used during classifier initialization and train-test splits
        :return: estimated value of c: probability that a positive example is labeled: c=p(s=1|y=1)
        """
        results = []
        X, y = self.create_training_data(P=P, U=U)
        # sampling k folds in a stratified manner
        skf = StratifiedKFold(n_splits=n_splits, random_state=random_state, shuffle=True)
        # get indices of train-test splits
        skf.get_n_splits(X, y)
        # iterating through indices
        for train_index, test_index in skf.split(X, y):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            model = self._choose_model(model=model, random_state=random_state)
            c = self._train_classifier_to_estimate_c(X_train, y_train, X_test, y_test, model)
            results.append(c)
        results = np.array(results)
        c_avg = np.average(results)
        c_std = np.std(results)
        print("Average of c: {:.2f}% ({:.2f})%".format(c_avg * 100, c_std * 100))
        self.c = c_avg
        return self.c

    def estimate_positive_count(self, P: list, U: list, n_splits=10, random_state=None):
        """
        Estimate the number of positive samples based on the C constant and number of positives. If
        :return: estimated number of positive samples from the unlabeled dataset
        """
        if self.c == 0:
            self.estimate_c(P=P, U=U, n_splits=n_splits, random_state=random_state)
        return int(len(list(P)) / self.c)
