from adze_modeler.model import BaseModel
import functools

class Simulation:
    app_name = 'adze project'

    def __init__(self, model:BaseModel):

        # the model creation class will be stored here
        self.model = model

        # These variables are reserved for the communication with the server.
        self._input = {}
        self._output = {}

        # These variables store the sections of the input data.
        self.cfg_simulation= {}
        self.cfg_model = {}
        self.cfg_tolerances = {}
        self.cfg_misc = {}

        # This dictionary stores the functions for the different simulations.
        self.simulations = {}

    def set_model(self, model:BaseModel):
        """
        Set a model class for the simulations.

        Parameters:
            model: This should be a subclass of the ModelBase class.
        """
        # TODO: check if model is class and not an object.

        assert issubclass(model, BaseModel), "model is not a BaseModel subclass."
        self.model = model

    def update_input(self):
        self.cfg_simulation = self._input['simulation']
        self.cfg_model = self._input['model']
        self.cfg_tolerances = self._input['tolerances']
        self.cfg_misc = self._input['misc']

        if self.cfg_simulation['type'] not in self.simulations.keys():
            raise ValueError(f"There is no simulation called {self.cfg_simulation['type']!r}")


        if 'exportname' in self.cfg_misc.keys():
            self.cfg_model['exportname'] = self.cfg_misc.pop('exportname')

        # TODO: decide model params.json
        # if self.cfg_tolerances['parameters']:
        #     for param_i in self.cfg_tolerances['parameters']:
        #         if param_i not in self.cfg_model.keys():
        #             raise ValueError(f'The model parameter {param_i!r} does not exist.')



    def run(self):
        """
        This is the main handler of an API call. After the input validation
        this function will execute the selected simulation and puts the results
        into the _output dictionary.
        """

        # get the simulation type from the simulation section.
        sim_type = self.cfg_simulation.pop('type')

        # If any parameter is present in the tolerances section, then a
        # tolerance analysis will be executed. Otherwise call the registered
        # function with the input arguments.
        if self.cfg_tolerances['parameters']:
           self.tolerance_analysis()
        else:
            self._output['res'] = self.simulations[sim_type](self.model,
                    self.cfg_model, self.cfg_simulation, self.cfg_misc)

    def tolerance_analysis(self):
        ...

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
            functools.wraps(func)
            def _wrapper(*arg, **kw):
                return func(*arg, **kw)
            return _wrapper

        return _decorator

sim = Simulation()
