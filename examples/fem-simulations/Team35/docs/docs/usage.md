# Usage

You can interact with the API and the digital twin behind it via simple make calls through code or via the `/apidocs` endpoint.

## Api calls by `/apidocs` endpoint

Calling the `/process_sim` endpoints, the interaction happens in a JSON file that is passed along with the call.

> TODO
> In the *fem-simulations* example we examine the Team35 benchmark problem where you can try simulations.
> There are more simulation types, but we use the **Basic** one to presenting the problem.

#### Input JSON File format

The input JSON file has several sections, these are the **simulation, model, tolerances, misc**. The only necessary section is the simulation, precisely the type field in the simulation section. All the other sections are optional.

```json
{
  "simulation": {
    "type": "default"
  },
  "model": {},
  "tolerances": {
    "type": "ff",
    "parameters": {},
    "variables": []
  },
  "misc": {
    "processes": 4,
    "cleanup": true
  }
}
```

In the `simulation` section you can specify the selected simulation as well as its additional parameters if it has any.

In the `model` section you can overwrite the provided model parameters.

On `tolerance` section you can make a tolerance analisys.

In the `misc` section, other variables are listed that are not tightly coupled with the digital twin. For example is a simulation is running in parallel you can specify the number of processes here.

### Use Basic simulation

Run a calculation with the given `model?` parameters. This simulation will give back
the [f1 function](#link-to-explain-f1-formula). If no input is given the simulation will use the default values.

The simulations and them default values defined in the **simulation.json** file.

The parameters of the simulation:

* **type**: simulation type
* **B0**: magnetic flux density in origo
* **x**: input param's vector *(consist of 10 elements)*

> TODO: define params meaning, create reference link

#### Api call with specifying parameters

> Note, that the 'x' vector require 10 elements!

Example input

```json
{
  "simulation": {
    "type": "default",
    "x": [7, 8, 9, 10, 11, 12, 13, 14, 15, 20],
    "B0": 3e-2
  }
}
```


#### Api call with pre-defined params, without specifying our parameters

Example input

```json
{
  "simulation": {
    "type": "basic"
  }
}
```

Example output

```json
{
  "res": {
    "f1": 0.0007730223408010202
  }
}
```


## Make calls through code
>TODO

## Simulation's platform

You can specify the model's platform.

There are several platforms that the application can manage.

The chosen simulation can be made by the selected platform.
The user can change the platform by rewrite a simple line in a code, and can export a simulation into Agros2D what previously was created by FEMM.

You can just give the suitable value to the `platform` variable:

* platform_femm
* platform_agros

```python
from digital_twin_distiller import (BaseModel,Snapshot,...)

class DistributedWinding(BaseModel):
    ...
    def setup_solver(self):
        ...
        platform = platform_femm
        # platform = platform_agros
        self.snapshot = Snapshot(platform)
...
```

> TODO: connection with model and simulation over the example, **model.json**

> TODO: running this **model.py** and the usage of the **/process_sim** endpoint

szimuláció pedig elkészíthető a kiválasztott programmal. Azaz a felhasználó egyetlen parancs módosításával képes például egy FEMM-ben megvalósított mágneses szimulációt Agros2D-be kiexportálni.

> TODO: FEMM, Agros2D mentions and reference link to them
