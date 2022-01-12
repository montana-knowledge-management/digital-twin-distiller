# 2. Usage

You can interact with the API and the digital twin behind it via simple make calls through code or via the `/apidocs`
endpoint.

## 2.1 Api calls by `/apidocs` endpoint

Calling the `/process_sim` endpoints, the interaction happens in a JSON file that is passed along with the call.

> TODO
> In the *fem-simulations* example we examine the Team35 benchmark problem where you can try simulations.
> There are more simulation types, but we use the **Basic** one to presenting the problem.

### 2.1.1 Input JSON File format

The input JSON file has several sections, these are the **simulation, model, tolerances, misc**. The only necessary
section is the simulation, precisely the `type` field in the simulation section. All the other sections are optional.

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

#### 2.1.1.1 Simulation section

You can specify the selected simulation as well as its additional parameters if it has any. The simulations and them
default values are defined in the `Team35/defaults/simulation.json` file.

``` json title="simulation.json"
{
  "default": {
      "B0": 2e-3,
      "x": [6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
  }
}
```

#### 2.1.1.2 Model section

You can overwrite the provided model parameters. The default model values can be found into
the `Team35/defaults/model.json` file.

``` json title="model.json"
{
  "x0": 1.0,
  "mw": 5
}
```

> TODO: what is $x_0$ and 'mw'

#### 2.1.1.3 Tolerance section

On this section you can make a tolerance analysis.
> TODO

#### 2.1.1.4 Misc section

In this section other variables are listed that are not tightly coupled with the digital twin. For example if a
simulation is running in parallel you can specify the number of processes here. The default mics values can be found
into the `Team35\defaults\mics.json` file.

```json title="mics.json"
{
  "processes": 4,
  "cleanup": true
}
```

### 2.1.2 Use Default simulation

Run a calculation with the given model parameters. This simulation will give back
the [f1 function](index.md#12-the-f_1-function). If no input is given the simulation will use the default values.

The simulations and them default values defined in the `simulation.json` file.

The parameters of the simulation:

| Parameter            | Value                               |
|:---------------------|:------------------------------------|
| type                 | simulation type                     |
| [B0](index.md#b0)    | magnetic flux density in the pole   |
| [x](index.md#radii)  | vector, consist of 10 radii params  |

#### 2.1.2.1 Api call with specifying parameters

This example has already demonstrated in the [main page](index.md#13-basic-usage-example-over-api-call)!

#### 2.1.2.2 Api call with pre-defined params, without specifying our parameters

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

## 2.2 Make calls through code

> TODO

## 2.3 Simulation's platform

The next step is to set up a FEM solver. Since there will be no exotic boundary conditions you can use either Agros2D or
FEMM. In this example we provide a setup for both of the software, and you can choose between them by changing the
highlighted line.

You can find more information about the **solver setup** on the [... page](#link-to-long-code-example)

You can set up a FEM solver.

There are several platforms that the application can manage.

The chosen simulation can be made by the selected platform. The user can change the platform by rewrite a simple line in
a code, and can export a simulation into Agros2D what previously was created by FEMM.

You can just give the suitable value to the `platform` variable:

* platform_femm
* platform_agros

```python
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

szimuláció pedig elkészíthető a kiválasztott programmal. Azaz a felhasználó egyetlen parancs módosításával képes például
egy FEMM-ben megvalósított mágneses szimulációt Agros2D-be kiexportálni.

> TODO: FEMM, Agros2D mentions and reference link to them
