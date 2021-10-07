from model import BaseModel
from abc import ABCMeta
import functools



class Simulation:
    app_name = 'adze project'
    simulations = {}

    def __init__(self, model=BaseModel):
        self.model = model

        # input dictionary contains the optimized parameters, tolerances and the setup strings
        self._input = {}
        self._output = {}

        self.cfg_simulation= {}
        self.cfg_model = {}

        self.simulations = {}

    def set_model(self, model):
        self.model = model

    def run(self):
        self._output['res'] = self.simulations['SIMPLE'](self.model, self.cfg_model)


    def register(self, name):
        def _decorator(func):
            self.simulations[name] = func
            functools.wraps(func)
            def _wrapper(*arg, **kw):
                return func(*arg, **kw)
            return _wrapper

        return _decorator

sim = Simulation()