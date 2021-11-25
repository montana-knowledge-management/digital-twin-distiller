import os
from abc import abstractmethod
from http.client import responses
from json import load

import tika.parser


class Reader:
    root_dir = os.path.dirname(os.path.abspath(os.path.join(__file__, "..")))

    @abstractmethod
    def read(self, file_name):
        pass


class JsonReader(Reader):
    def __init__(self):
        super().__init__()

    def read(self, file_name):
        content = None

        try:
            with open(file_name, encoding="utf-8") as json_file:
                content = load(json_file)
        except ValueError:
            # handling unexpected bom headers
            with open(file_name, encoding="utf-8-sig") as json_file:
                content = load(json_file)

        except OSError:
            pass

        return content


class PdfReader(Reader):
    """
    This reader gives back a dictionary, which contains the text content under the standard Text key, then the
    code gives back a status key and the other keys which found in the dictionary.

    :param file_name:
    :return:
    """

    ok = responses[200]
    error = responses[503]
    text_key = "Text"
    status = "status"
    meta = "metadata"

    def __init__(self):
        super().__init__()
        self.server_init()

    def server_init(self):
        tika.initVM()

    def read(self, file_name: str = None) -> dict:
        ret_dict = {self.status: self.error}

        # the tika can't work with posixpath, we are forcing to converse these path into strings
        file_name = str(file_name)
        try:
            parsed = tika.parser.from_file(file_name)
            if parsed.get(self.meta):
                ret_dict.update(parsed.get(self.meta))
            ret_dict[self.text_key] = parsed.get("content")
            ret_dict[self.status] = responses[parsed.get(self.status)]

        except OSError:
            pass

        return ret_dict


class TextReader(Reader):
    """
    The class reads a given text this function should be invoked via the
    document opened function of the Abstract Project.
    """

    def read(self, file_name):
        with open(file_name, encoding="utf-8") as text_file:
            content = text_file.read()
            return content
