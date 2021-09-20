# Benchmark TEAM Problem for coil optimization

> Di Barba, P.; Mognaschi, M.E.; Lowther, D.A.; Sykulski, J.K.  A benchmark TEAM
> problem for multi-objective Pareto optimization of electromagnetic
> devices. IEEE Transactions on Magnetics, 2018, 54, 1–4.

In this application, we will optimize one of the examples with the help of the Artap framework,
namely the [distributed winding coil](../../examples/DistributedWinding/readme.md) example.

## Symmetric case

```python
class CoilOptimizationProblem(Problem):
    def set(self):
        self.name = "Symmetric Distributed Winding Test Problem"

        self.parameters = [
            {"name": "r0", "bounds": [5.5, 20]},
            {"name": "r1", "bounds": [5.5, 20]},
            {"name": "r2", "bounds": [5.5, 20]},
            {"name": "r3", "bounds": [5.5, 20]},
            {"name": "r4", "bounds": [5.5, 20]},
            {"name": "r5", "bounds": [5.5, 20]},
            {"name": "r6", "bounds": [5.5, 20]},
            {"name": "r7", "bounds": [5.5, 20]},
            {"name": "r8", "bounds": [5.5, 20]},
            {"name": "r9", "bounds": [5.5, 20]},
        ]

        self.costs = [
            {"name": "f_1", "criteria": "minimize"},
            {"name": "f_2", "criteria": "minimize"},
            {"name": "f_3", "criteria": "minimize"},
        ]

    def evaluate(self, individual):
        perturbation = 0.5  # mm
        X = individual.vector
        Xp = [xi + perturbation for xi in X]
        Xn = [xi - perturbation for xi in X]
        models = [SymmetircModel(X), SymmetircModel(Xn), SymmetircModel(Xp)]

        F1 = math.inf
        F2 = math.inf
        F3 = math.inf
        try:

            with multiprocessing.Pool(processes=3) as pool:
                res, resn, resp = pool.map(execute_model, models)

            Bz = [pointvalue[2] for pointvalue in res["Bz"]]  # [x, y, Bz(x, y)]
            Br = [pointvalue[2] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
            xi = [pointvalue[0] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
            yi = [pointvalue[1] for pointvalue in res["Br"]]  # [x, y, Br(x, y)]
            nb_nodes = res["nodes"]

            # Calculate F1
            B0 = 2e-3
            F1 = max(map(lambda Bz_i: abs(Bz_i - B0), Bz))

            # Calculate F2

            Bzp = [pointvalue[2] for pointvalue in resp["Bz"]]
            Brp = [pointvalue[2] for pointvalue in resp["Br"]]

            Bzn = [pointvalue[2] for pointvalue in resn["Bz"]]
            Brn = [pointvalue[2] for pointvalue in resn["Br"]]

            deltaBpz = map(operator.abs, map(operator.sub, Bzp, Bz))
            deltaBpr = map(operator.abs, map(operator.sub, Brp, Br))
            deltaBp = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBpz, deltaBpr))

            deltaBnz = map(operator.abs, map(operator.sub, Bzn, Bz))
            deltaBnr = map(operator.abs, map(operator.sub, Brn, Br))
            deltaBn = map(math.sqrt, map(lambda a, b: a ** 2 + b ** 2, deltaBnz, deltaBnr))

            F2 = max(map(operator.add, deltaBp, deltaBn))

            # Calculate F3
            F3 = sum(X)

            print(F1, F2, F3)
            return [F1, F2, F3]
        except Exception as e:
            return [math.inf, math.inf, math.inf]
```
