from abc import abstractmethod
from adze_modeler.modelpiece import ModelPiece
import string
from random import choices
from shutil import rmtree
from adze_modeler.snapshot import Snapshot
from uuid import uuid4

class Model:
    def __init__(self):
        self.ingredients = dict()
        self.snapshot: Snapshot = None
        self.model_id = str(uuid4())


    def add_piece(self, p: ModelPiece):
        self.ingredients[p.name] = p

    @abstractmethod
    def build(self, p: dict):
        ...


    def cleanup(self):
        try:
            rmtree(self.model_path)

        except Exception as e:
            print(e)

