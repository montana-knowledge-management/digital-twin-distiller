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

    def update_input(self):
        self.cfg_simulation = self._input['simulation']
        self.cfg_model = self._input['model']
        self.cfg_misc = self._input['misc']

        if self.cfg_simulation['type'] not in self.simulations.keys():
            raise ValueError(f"There is no simulation called {self.cfg_simulation['type']!r}")



    def run(self):
        sim_type = self.cfg_simulation.pop('type')
        self._output['res'] = self.simulations[sim_type](self.model, self.cfg_model, self.cfg_simulation)


    def register(self, name):
        def _decorator(func):
            self.simulations[name] = func
            functools.wraps(func)
            def _wrapper(*arg, **kw):
                return func(*arg, **kw)
            return _wrapper

        return _decorator

sim = Simulation()
