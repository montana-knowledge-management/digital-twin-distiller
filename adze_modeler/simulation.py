from model import BaseModel
from abc import ABCMeta
from abc import abstractmethod


class Simulation(metaclass=ABCMeta):
    app_name = 'adze project'

    def __init__(self, model=BaseModel):
        self.model = model

        # input dictionary contains the optimized parameters, tolerances and the setup strings
        self._input = {}
        self._output = {}

    def set_model(self, model):
        self.model = model

    def run(self):
        ...
