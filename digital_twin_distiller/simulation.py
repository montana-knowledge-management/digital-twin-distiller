import functools
import json
import operator as op
import pprint
from typing import Dict, Sequence

from digital_twin_distiller.doe import *
from digital_twin_distiller.model import BaseModel
from digital_twin_distiller.modelpaths import ModelDir


class Simulation:
    app_name = "adze project"

    def __init__(self, model: BaseModel = None):

        # the model creation class will be stored here
        self.model = model

        # These variables are reserved for the communication with the server.
        self._input = {}
        self._output = {}

        # These variables store the sections of the input data.
        self.cfg_simulation = {}
        self.cfg_model = {}
        self.cfg_tolerances = {}
        self.cfg_misc = {}

        # This dictionary stores the functions for the different simulations.
        self.simulations = {}

    def set_model(self, model: BaseModel):
        """
        Set a model class for the simulations.

        Parameters:
            model: This should be a subclass of the ModelBase class.
        """
        # TODO: check if model is class and not an object.

        assert issubclass(model, BaseModel), "model is not a BaseModel subclass."
        self.model = model

    def _load_defaults(self):

        sim_type = self._input["simulation"]["type"]

        file_sim = ModelDir.DEFAULTS / "simulation.json"
        assert file_sim.exists(), f"Default simulation.json does not exist @ {file_sim.resolve()}"
        with open(file_sim) as f:
            default_cfg = dict(json.load(f))
            if sim_type not in default_cfg.keys():
                raise ValueError(f"There is no simulation called {self.cfg_simulation['type']!r}")

            self.cfg_simulation = default_cfg[sim_type]

        file_model = ModelDir.DEFAULTS / "model.json"
        assert file_model.exists(), "Default model.json does not exist."
        with open(file_model) as f:
            self.cfg_model = dict(json.load(f))

        file_misc = ModelDir.DEFAULTS / "misc.json"
        assert file_misc.exists(), "Default misc.json does not exist."
        with open(file_misc) as f:
            self.cfg_misc = dict(json.load(f))

    def update_input(self):
        self._load_defaults()

        # overwrite the default configs with the input
        self.cfg_simulation.update(self._input["simulation"])
        self.cfg_model.update(self._input["model"])
        self.cfg_tolerances.update(self._input["tolerances"])
        self.cfg_misc.update(self._input["misc"])

        if "exportname" in self.cfg_misc.keys():
            self.cfg_model["exportname"] = self.cfg_misc.pop("exportname")

        if self.cfg_tolerances["parameters"]:
            for param_i in self.cfg_tolerances["parameters"]:
                if param_i not in self.cfg_model.keys():
                    raise ValueError(f"The model parameter {param_i!r} does not exist.")

    def run(self):
        """
        This is the main handler of an API call. After the input validation
        this function will execute the selected simulation and puts the results
        into the _output dictionary.
        """

        # get the simulation type from the simulation section.
        sim_type = self.cfg_simulation["type"]

        # If any parameter is present in the tolerances section, then a
        # tolerance analysis will be executed. Otherwise call the registered
        # function with the input arguments.
        if self.cfg_tolerances["parameters"]:
            self.tolerance_analysis()
        else:
            self._output["res"] = self.simulations[sim_type](
                self.model, self.cfg_model, self.cfg_simulation, self.cfg_misc
            )

    def tolerance_analysis(self):
        """
        Execute a tolerance analysis with the selected parameters on a simulation.
        """
        sim_type = self.cfg_simulation["type"]
        doe_type = self.cfg_tolerances["type"]
        variables = self.cfg_tolerances["variables"]
        parameter_names = tuple(self.cfg_tolerances["parameters"].keys())
        parameter_tolerances = tuple(self.cfg_tolerances["parameters"].values())
        original_values = tuple(self.cfg_model[pi] for pi in parameter_names)

        result = self.simulations[sim_type](self.model, self.cfg_model, self.cfg_simulation, self.cfg_misc)
        self._output["res"] = result
        self._output["tolerances"] = {}
        yref = self._format_result(result)

        # get the designs
        designs = []
        if doe_type == "ff":
            designs = doe_fullfact([3] * len(parameter_names))
        elif doe_type == "pb":
            designs = doe_pbdesign(len(parameter_names))
        elif doe_type == "bb":
            designs = doe_bbdesign(len(parameter_names))
        elif doe_type == "ccf":
            designs = doe_ccf(len(parameter_names))

        Y = []
        for di in designs:
            dX = map(op.mul, di, parameter_tolerances)
            dX = map(op.add, dX, original_values)
            X = dict(zip(parameter_names, dX))
            self.cfg_model.update(X)

            results = self.simulations[sim_type](self.model, self.cfg_model, self.cfg_simulation, self.cfg_misc)

            Y.append(self._format_result(results))

        # print(pprint.pformat(Y, indent=2, compact=True))

        for var_i in variables:
            # self._output['tolerances'][var_i] = {'upper': None, 'lower': None, 'upper_res': None, 'lower_res': None}
            self._output["tolerances"][var_i] = {"upper": None, "lower": None}
            for yi in Y:
                delta = list(map(op.sub, yi[var_i], yref[var_i]))
                delta.sort(key=lambda di: abs(di))
                yi["delta"] = delta[-1]

            Y.sort(key=op.itemgetter("delta"))
            self._output["tolerances"][var_i]["upper"] = Y[-1]["delta"]
            # self._output['tolerances'][var_i]['upper_res'] = Y[-1]
            self._output["tolerances"][var_i]["lower"] = Y[0]["delta"]
            # self._output['tolerances'][var_i]['lower_res'] = Y[0]

    def _format_result(self, results) -> dict:
        """
        Format any result into a standard dictionary format. In this format all dictionary keys have a list as a value.

        Examples:
            - results = {'a':x1, 'b':[c1,c2,c3]} -> {'a': [x1], 'b':[c3]}
                - note: point coordinates will be discarded (c1, c2)

            - results = [{'a': x1, 'b': [c11, c12, c13]},
                         {'a': x2, 'b': [c21, c22, c23]},
                         {'a': x3, 'b': [c31, c32, c33]}]
                         ->
                         {'a':[x1, x2, x3], 'b':[c13, c23, c33]}

        Parameters:
            results: It can be a single dictionary or a list of dictionaries.
        """
        # this variable will be returned
        y = {}
        # check if the input is a list
        if not isinstance(results, Sequence):
            # if not a list transform it
            results = [results]

        # get the dictionary keys
        # NOTE: Assume all dictionaries have the same set of keys.
        variables = results[0].keys()
        for key in variables:
            # extract the values bounded to key from all the dictionaries
            y[key] = list(map(op.itemgetter(key), results))

        # If any keys have a value of list of list, then change it to list of
        # numbers by extracting the last element from the sublists.
        # NOTE: this is the place where the coordinates of the point values are
        #       getting discarded
        for key, value in y.items():
            if isinstance(value[0], Sequence):
                y[key] = list(map(op.itemgetter(-1), value))

        return y

    def register(self, name):
        """
        Register an ordinary function as a simulation. The function should have
        the following signature: function_name(model, modelparams, simprams, miscparams)
        and should return a dict with the results or a list of dicts.

        Parameters:
            name: The name of the simulation. This name will be used in the
                  json API call to identify the simulation.
        """

        # Because register is a parametric decorator, an inside function should
        # be defined. This function will be returned as the real decorator.
        def _decorator(func):
            # register the function in the simulations dictionary
            self.simulations[name] = func

            # this is the wrapper function around the original.  This section
            # is needed to use the input function as a normal function.
            @functools.wraps(func)
            def _wrapper(*arg, **kw):
                return func(*arg, **kw)

            return _wrapper

        return _decorator


sim = Simulation()
