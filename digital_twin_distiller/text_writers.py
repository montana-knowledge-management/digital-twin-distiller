from abc import abstractmethod
from json import dump


class Writer:
    def __init__(self):
        ...

    @abstractmethod
    def write(self, data_dictionary, file_name):
        pass


class JsonWriter(Writer):
    def write(self, data_dictionary, file_name, encoding="utf-8"):
        success = False
        """ Writing out the extracted data to a json file """
        try:
            with open(file_name, "w", encoding=encoding) as outfile:
                dump(data_dictionary, outfile, indent=4, ensure_ascii=False)
                success = True
        except OSError:
            pass
        return success


class TextWriter(Writer):
    def write(self, text: str, file_name):
        """Writing out text into file"""
        success = False
        try:
            with open(file_name, "w") as text_file:
                text_file.write(f"{text}")
                success = True
        except OSError:
            pass
        return success
