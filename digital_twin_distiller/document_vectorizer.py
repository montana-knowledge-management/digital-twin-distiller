import numpy
import numpy as np
from digital_twin_distiller.ml_project import AbstractTask
import fasttext
from gensim.models import FastText, Word2Vec
from importlib_resources import files
from digital_twin_distiller.text_readers import JsonReader
from digital_twin_distiller.text_writers import JsonWriter
from sklearn.feature_extraction.text import TfidfVectorizer
from tqdm import tqdm
from gensim.models.doc2vec import TaggedDocument, Doc2Vec


class DocumentVectorizer(AbstractTask):
    def __init__(self):
        self.vocabulary = None
        self.gensim_model = None
        self.fasttext_model = None
        self.doc2vec_model = None
        self.vectors_dict = {}
        self.original_text = None
        self.augmented_text = []
        self.idf = []

    def load_gensim_model(self, model_path):
        """
        Loads gensim FastText or Word2Vec models.
        :param fasttext_model_path: path for pretrained model
        :return:
        """
        try:
            self.gensim_model = Word2Vec.load(model_path)
        except Exception:
            self.gensim_model = FastText.load(model_path)

    def load_fasttext_model(self, fasttext_model_path):
        """
        Loads models from https://fasttext.cc/docs/en/crawl-vectors.html. Only bin format is supported!
        :param fasttext_model_path: path for pretrained model in bin format is supported.
        :return:
        """
        self.fasttext_model = fasttext.load_model(fasttext_model_path)

    def load_vectors_dictionary(self, path_to_dict):
        """
        Load previously created vectors dictionary containing vectors for words from json file.
        :param path_to_dict:
        :return:
        """
        reader = JsonReader()
        loaded_json = reader.read(path_to_dict)
        casted_vectors = {}
        # casting from list to numpy array
        for word, vector in loaded_json.items():
            casted_vectors[word] = np.array(vector)

        self.vectors_dict = casted_vectors

    def save_vectors_dictionary(self, path_to_dict):
        """
        Writes the content of self.vectors_dict to json file.
        :param path_to_dict: path to save
        :return: None
        """
        writer = JsonWriter()
        casted_vocabulary = {}
        # casting from numpy array to list
        for key, value in self.vectors_dict.items():
            casted_vocabulary[key] = value.tolist()
        writer.write(casted_vocabulary, path_to_dict)

    def build_vocab(self, text, **kwargs):
        """
        For faster processing, a vocabulary has to be created to perform similarity actions only once per token and not
        multiple times per document.
        :param text: list of strings, non-tokenized
        :param kwargs: parameters of a CountVectorizer object e.g. lowercase
        :return: None
        """
        if not kwargs:
            kwargs = {"lowercase": False}
        count_vect = TfidfVectorizer(**kwargs)
        count_vect.fit_transform(text)
        self.vocabulary = count_vect.vocabulary_
        self.idf = count_vect.idf_

    def build_vectors_dict(self, mode="fasttext"):
        """
        Must be called when the vocabulary has been built.
        :return: A dictionary containing most similar finds. Keys are the members of the vocabulary,
        """
        if mode == "gensim":
            if not self.gensim_model:
                raise ValueError("Missing loaded gensim model. Please load gensim model!")
            known_words = set(self.vectors_dict.keys())
            for word in tqdm(self.vocabulary):
                if not {word}.intersection(known_words):
                    try:
                        word_vector = self.gensim_model.wv[word]
                        self.vectors_dict[word] = word_vector
                    except KeyError:
                        print("The word {} was missing from the gensim model, not putting into vectors dict.".format(
                            word))

        elif mode == "fasttext":
            if not self.fasttext_model:
                raise ValueError("Missing loaded fasttext model. Please load fasttext model!")
            known_words = set(self.vectors_dict.keys())
            for word in tqdm(self.vocabulary):
                # using sets improves execution speed
                if not {word}.intersection(known_words):
                    word_vector = self.fasttext_model.get_word_vector(word)
                    self.vectors_dict[word] = word_vector
        else:
            raise ValueError("Wrong mode given! Please choose from 'gensim' or 'fasttext'!")

    @staticmethod
    def cosine_similarity(vector_1: numpy.ndarray, vector_2: numpy.ndarray):
        cos_sim = np.dot(vector_1, vector_2) / (np.linalg.norm(vector_1) * np.linalg.norm(vector_2))
        return cos_sim

    def keep_n_most_similar_to_average(self, tokenized_document: list, nr_of_words_to_keep=10, mode="average"):
        """
        Keeps the words that are most similar to the document vector calculated by word vector average or idf weighted
        average.
        :param mode: average: document vetor = average of vectors in text
                     idf_weighted: document vetor = idf weighted average of vectors in text
        :param tokenized_document: list of tokens
        :param nr_of_words_to_keep: int: number of words to keep that are the closest to the average document vector
        """
        if mode == "average":
            avg_vector, vectors = self.calculate_average(tokenized_document)
        elif mode == "idf_weighted":
            avg_vector, vectors = self.calculate_idf_weighted_average(tokenized_document)
        else:
            raise ValueError('The mode is currently not implemented, please choose from "average" and "idf_weighted"!')
        similarities = [(self.cosine_similarity(avg_vector, vector), vector) for vector in vectors]
        similarities = sorted(similarities, key=lambda x: x[0])[:nr_of_words_to_keep]
        similarities = [vec[1] for vec in similarities]
        similarities = np.array(similarities)
        return np.mean(similarities, axis=0)

    def calculate_average(self, tokenized_document: list):
        """
        Gets the vectors from the model and calculates the average of vectors for each document.
        :param tokenized_document: list of tokens of one document
        :returns: average of vectors, and the vectors themselves
        """
        vectors = [self.fasttext_model.get_word_vector(token) for token in tokenized_document if token]
        vectors = np.array(vectors)
        return np.mean(vectors, axis=0), vectors

    def calculate_idf_weighted_average(self, tokenized_document: list):
        """
        Gets the vectors from the model and calculates the idf weighted average of vectors for each document.
        If the token cannot be found in the vocabulary, the factor remains 1.0.
        :param tokenized_document: list of tokens of one document
        :returns: idf weighted average of vectors, and the vectors themselves
        """
        vectors = []
        for token in tokenized_document:
            vector = self.fasttext_model.get_word_vector(token)
            if self.vocabulary.get(token) is not None:
                factor = self.idf[self.vocabulary.get(token)]
            else:
                factor = 1.0
            vector = vector * factor
            vectors.append(vector)
        # ids = [self.idf[self.vocabulary.get(token)] if self.vocabulary.get(token) else 1.0 for token in tokenized_document]
        # print(ids)
        # vectors = [self.fasttext_model.get_word_vector(token) for token in tokenized_document if token]
        vectors = np.array(vectors)
        return np.mean(vectors, axis=0), vectors

    def calculate_doc2vec(self, tokenized_document: list):
        """
        Gets the vectors from the doc2vec model.
        :param tokenized_document: list of tokens of one document
        :returns: doc2vec vector
        """
        # checking whether trained model is available
        if not self.doc2vec_model:
            raise ValueError("Missing Doc2Vec model!")
        return self.doc2vec_model.infer_vector(tokenized_document)

    def train_doc2vec_model(self, corpus, model_path_to_save=None, **kwargs):
        """
        Trains doc2vec model on a given corpus and saves it as long as the path to save is specified.
        :param corpus: list of documents in  list of strings format
        :param model_path_to_save: default: None, if given, the trained model will be saved automatically to this path
        :param kwargs: arguments that are the same as it would be for a doc2vec model for full set of parameters please
                       visit: https://radimrehurek.com/gensim/models/doc2vec.html
            examples for parameters:
            :param vector_dim: required dimension of the doc2vec model
            :param epochs: epochs to be used during training
            :param min_count: mininum required count of tokens
        """
        train_corpus = []
        for idx, document in enumerate(corpus):
            train_corpus.append(TaggedDocument(document.split(), [idx]))
        model = Doc2Vec(**kwargs)
        model.build_vocab(train_corpus)
        model.train(train_corpus, total_examples=model.corpus_count, epochs=model.epochs)
        self.doc2vec_model = model
        if model_path_to_save:
            model.save(model_path_to_save)
        return model

    def load_doc2vec_model(self, path_to_doc2vec_model):
        """
        Loads pretrained Doc2Vec model.
        """
        self.doc2vec_model = Doc2Vec.load(path_to_doc2vec_model)
        return self.doc2vec_model

    def run(self, tokenized_document: list, mode="average"):
        """
        Performs the required vectorization on the given documents.
        :param tokenized_document: list of list of tokens of a document
        :param mode: available modes:
                                "average": returns average vector for the documents
                                "idf_weighted": returns idf weighted average vector for the documents
                                "doc2vec": returns the doc2vec vector for the documents
        :returns: list of vectors for the documents
        """
        if isinstance(tokenized_document, list) and isinstance(tokenized_document[0], str):
            tokenized_document = [tokenized_document]
        return_vectors = []
        for document in tokenized_document:
            if mode == "average":
                return_vectors.append(self.calculate_average(document)[0])
            elif mode == "idf_weighted":
                return_vectors.append(self.calculate_idf_weighted_average(document)[0])
            elif mode == "doc2vec":
                return_vectors.append(self.calculate_doc2vec(document))
        return return_vectors
