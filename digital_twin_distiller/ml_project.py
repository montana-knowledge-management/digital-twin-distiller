import logging
import os
import socket
from abc import abstractmethod

import pkg_resources
import requests
from importlib_resources import files

import digital_twin_distiller.text_readers as rdr
from digital_twin_distiller.keywords import JSON, PDF, TXT
from digital_twin_distiller.text_writers import JsonWriter

supported_extensions = [JSON, TXT, PDF]


class MachineLearningProject:
    log = logging.getLogger(__name__)
    # This folder contains the index page for the api endpoint
    # template_folder = files("digital_twin_distiller") / "resources" / "templates"
    # static_folder = files("digital_twin_distiller") / "resources" / "static"

    port = 9099
    debug = False
    cached_subtasks = {}

    def __init__(self, app_name="Distiller", no_cache=False):
        self._output_data = (
            []
        )  # every data file is represented by a dictionary, a bulk input can be imported as a list of dicts
        self._temp_data = []
        self._input_data = []
        self.options = {"res_dir": pkg_resources.resource_dir}  # the ini file options list can be stored here
        self.pipe_document_extractors = []
        self.app_name = app_name
        self.host = self.get_ip()[0]

        if not no_cache:
            self.cache()

    @abstractmethod
    def custom_input(self):
        pass

    def add_single_input(self, json_data: dict):
        self._input_data = [json_data]

    def get_single_output(self):
        if len(self._output_data):
            return self._output_data[0]

    def get_ip(self):
        # helper to find own ip address and opens a port near thats
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
        s.connect(("<broadcast>", 0))
        ip = s.getsockname()[0]
        s.close()

        # add the port number
        host = ip
        ip = f"{ip}:{self.port}"

        return (host, ip)

    @abstractmethod
    def cache(self):
        pass

    # TODO: this part is going to be replaced by the S3 bucket solution.
    # def http_input(self, http_input, file_format, reader=None, output_directory="/tmp/"):
    #     """
    #     This function should be finished and discussed! Not finished yet!
    #     """
    #
    #     response = requests.get(http_input)
    #     filename = response.url.split("/")[-1]
    #
    #     if file_format == PDF:
    #         path = os.path.join(output_directory, filename)
    #
    #         with open(path, "wb") as file:
    #             file.write(response.content)
    #         if reader:
    #             data = reader().read(path)
    #             JsonWriter().write(data, "tmp.json")

    def bulk_input_directory(self, directory_name, extension=JSON, key=None):
        """
        Reads the data into a dictionary, where the key is the filename than the
        """

        if extension not in supported_extensions:
            raise ValueError(f"{extension} - extension is not supported.")

        for file in os.listdir(directory_name):
            inp_dict = None
            if file.endswith(extension):
                file_path = os.path.join(directory_name, file)
                if extension == JSON:
                    inp_dict = self.open_data_file(rdr.JsonReader(), file_path, key)

                if extension == PDF:
                    inp_dict = self.open_data_file(rdr.PdfReader(), file_path, key)

                if extension == TXT:
                    inp_dict = self.open_data_file(rdr.TextReader(), file_path, key)

                if inp_dict:
                    self._input_data.append(inp_dict)

    def open_data_file(self, reader, file_name, key=None):
        """
        Stores the input data in  a dictionary, or the content under a dictionary key.

        Possible usage of the function in different example cases:

        - a text file: reads the text under the defined key -> {key: text file}
        - json: reads the json into a dictonary
        - json under a key, reads the
        - xml : reads the xml file into a dictionary
        """

        input_data = {}

        if not key:
            input_data = reader.read(file_name)

        else:
            input_data[key] = reader.read(file_name)

        return input_data

    def write_data_file(self, writer, output_file):
        writer.write(self._output_data, output_file)

    @abstractmethod
    def run(self):
        ...

    @abstractmethod
    def train(self):
        ...

    @abstractmethod
    def preprocess(self):
        ...


class AbstractTask:
    def __init__(self, project_inp=None, cache=None):
        self._input_data = project_inp
        self._output_data = {}  # task has its own output, which can be merged with other tasks output
        self._temp_data = {}  # task has its own temp
        self.sub_tasks = []  # list of the contained sub-tasks
        self.define_subtasks(cache=cache)

    @abstractmethod
    def define_subtasks(self, cache=None):
        ...

    def simple_execution(self):
        """Iterates over the defined"""
        ret_list = []
        for item in self._input_data:
            for s_t in self.sub_tasks:
                ret_list.append(s_t.run(item))
        return ret_list

    @abstractmethod
    def execute(self, input_data):
        ...

    def __call__(self, sub_tasks: list = None):
        # Simple and straightforward mechanism for executor when it executes the functions step by step
        self.sub_tasks = sub_tasks
        results = self.simple_execution()
        return results


class AbstractSubTask:
    @abstractmethod
    def run(self, input_data):
        ...

    def __call__(self, inp_dict):
        return self.run(inp_dict)


class Classifier(AbstractSubTask):
    def __init__(self):
        super().__init__()

    @abstractmethod
    def load_model(self):
        pass


class ExtractorAbstract(AbstractSubTask):
    def __init__(self):
        super().__init__()
        self.define_options()

    @abstractmethod
    def define_options(self):
        ...


class PreProcessorAbstract(AbstractSubTask):
    def __init__(self):
        super().__init__()
