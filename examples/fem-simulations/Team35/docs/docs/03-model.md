# Model creation

There are more than one way to create a model for a problem, but the simplest one is to just use the built in `BaseModel` class. This class implements most of the boilerplate code that is needed for the calculations. It sets up the proper paths and automatically builds, executes and retrieves the results.

**The users only required to fill the provided abstract methods to get a fully functional model.**

## Create a model using `BaseModel` class

Let's begin the model creation with subclassing the BaseModel class. To use this class we have to fill the {==highlighted==} methods. The class will have 1 input, which is the X list that holds the distance from the z axis for the different turns. In our geometry simulation a 10 element list is needed.

> !!! TODO: this came from SymmetircModel example, but we don't specify the symmetric and asymmetric cases
>
> Idea: refactor needed

``` python hl_lines="18 21 24 27 30"

class DistributedWinding(BaseModel):
    def __init__(self, X:list, exportname=None):
        # check for proper length
        assert len(X) == 10

        # Be aware of mutable arguments !
        # Use the copy of the argument to prevent bugs.
        self.X = X.copy()

        # call the base class constructor to set up paths
        super(SymmetircModel, self).__init__(exportname=exportname)

        # if you want to use different paths then your changes go here

        # create the directories
        self._init_directories()

    def setup_solver(self):
        ...

    def define_boundary_conditions(self):
        ...

    def define_materials(self):
        ...

    def add_postprocessing(self):
        ...

    def build_geometry(self):
        ...
```
