import numpy
import numpy as np
from digital_twin_distiller.ml_project import AbstractTask
import fasttext
from gensim.models import FastText
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
        self.vectors_dict = {}
        self.original_text = None
        self.augmented_text = []
        self.idf = []

    # def define_subtasks(self, cache=None):
    #     """
    #     Loads pretrained gensim language model (e.g. FastText, Word2Vec, etc.)
    #     :param cache: enables that loading is performed only once
    #     :return:
    #     """
    #     self.fasttext_model = fasttext.load_model(
    #         str(files("distiller") / "resources" / "augmentation" / "cc.hu.2.bin")
    #     )

    def load_gensim_model(self, model_path):
        """
        Loads gensim FastText models. Only bin format is supported!
        :param fasttext_model_path: path for pretrained model in bin format is supported.
        :return:
        """
        self.gensim_model = FastText.load(model_path)

    def load_fasttext_model(self, fasttext_model_path):
        """
        Loads models from https://fasttext.cc/docs/en/crawl-vectors.html. Only bin format is supported!
        :param fasttext_model_path: path for pretrained model in bin format is supported.
        :return:
        """
        self.fasttext_model = fasttext.load_model(fasttext_model_path)

    def load_most_similar_dictionary(self, path_to_dict):
        """
        Load previously created dictionary containing most similar tokens from json file.
        :param path_to_dict:
        :return:
        """
        reader = JsonReader()
        loaded_json = reader.read(path_to_dict)
        self.vectors_dict = loaded_json

    def save_most_similar_dictionary(self, path_to_dict):
        """
        Write most similar dictionary to json file.
        :param path_to_dict: path to save
        :return:
        """
        writer = JsonWriter()
        writer.write(self.vectors_dict, path_to_dict)

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
        # if mode == "gensim":
        #     if not self.gensim_model:
        #         raise ValueError("Missing loaded gensim model. Please load gensim model!")
        #     for word in self.vocabulary:
        #         if word in self.gensim_model.wv.key_to_index:
        #             if word not in self.most_similar_dict:
        #                 self.most_similar_dict[word] = self.gensim_model.wv.most_similar(
        #                     word, topn=self.topn_most_similar
        #                 )
        #         else:
        #             self.most_similar_dict[word] = word
        if mode == "fasttext":
            if not self.fasttext_model:
                raise ValueError("Missing loaded fasttext model. Please load fasttext model!")
            for word in tqdm(self.vocabulary):
                # using sets improves execution speed
                if not {word}.intersection(set(self.vectors_dict.keys())):
                    # if word not in self.most_similar_dict:
                    word_vector = self.fasttext_model.get_word_vector(word)

                    self.vectors_dict[word] = word_vector
        else:
            raise ValueError("Wrong mode given! Please choose from 'gensim' or 'fasttext'!")

    @staticmethod
    def cosine_similarity(vector_1: numpy.ndarray, vector_2: numpy.ndarray):
        cos_sim = np.dot(vector_1, vector_2) / (np.linalg.norm(vector_1) * np.linalg.norm(vector_2))
        return cos_sim

    def keep_n_most_similar_to_average(self, tokenized_document: list, nr_of_words_to_keep=10, mode="average"):
        if mode == "average":
            avg_vector, vectors = self.calculate_average(tokenized_document)
        elif mode == "idf_weighted":
            avg_vector, vectors = self.calculate_idf_weighted_average(tokenized_document)
        similarities = [(self.cosine_similarity(avg_vector, vector), vector) for vector in vectors]
        similarities = sorted(similarities, key=lambda x: x[0])[:nr_of_words_to_keep]
        similarities = [vec[1] for vec in similarities]
        similarities = np.array(similarities)
        return np.mean(similarities, axis=0)

    def calculate_average(self, tokenized_document: list):
        # vectors = []
        # for token in tokenized_document:
        #     if token:
        #         print(token, self.vectors_dict.get(token))
        #         vectors.append()
        vectors = [self.fasttext_model.get_word_vector(token) for token in tokenized_document if token]
        vectors = np.array(vectors)
        return np.mean(vectors, axis=0), vectors

    def calculate_idf_weighted_average(self, tokenized_document: list):
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
        # checking whether trained model is available
        if not self.doc2vec_model:
            raise ValueError("Missing Doc2Vec model!")
        return self.doc2vec_model.infer_vector(tokenized_document)

    def train_doc2vec_model(self, corpus, vector_dim, epochs, min_count, model_path_to_save=None, **kwargs):
        train_corpus = []
        for idx, document in enumerate(corpus):
            train_corpus.append(TaggedDocument(document.split(), [idx]))
        model = Doc2Vec(vector_size=vector_dim, min_count=min_count, epochs=epochs, **kwargs)
        model.build_vocab(train_corpus)
        model.train(train_corpus, total_examples=model.corpus_count, epochs=model.epochs)
        self.doc2vec_model = model
        if model_path_to_save:
            model.save(model_path_to_save)
        return model

    def load_doc2vec_model(self, path_to_doc2vec_model):
        self.doc2vec_model = Doc2Vec.load(path_to_doc2vec_model)

    def run(self, tokenized_document: list, mode="average"):
        if mode == "average":
            return self.calculate_average(tokenized_document)[0]
        elif mode == "idf_weighted":
            return self.calculate_idf_weighted_average(tokenized_document)[0]
        elif mode == "doc2vec":
            return self.calculate_doc2vec(tokenized_document)
