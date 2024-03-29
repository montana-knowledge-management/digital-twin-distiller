# 1. Distributed winding coil

Summary
-------

The proposed model can calculate one layout for the TEAM35 Test problem. You can select from the two integrated FEM
solvers to compute the results of the problem, a 10 turned solenoid, where the position of the different turns can be optimized to
generate a homogenous magnetic field distribution in the coil.

## 1.1 Basic usage example over api call

Endpoint: `/apidocs/process_sim`

In this example we use the default simulation without any additional [JSON section parameters](02-usage.md#211-input-json-file-format).
We specify only the `simulation` block.

The `x` vector input variable is represents a simple `radii` (turn distance from z axis), and it's size has to be 10.

### 1.1.1 Example input

```json
{
  "simulation": {
    "type": "default",
    "x": [7, 8, 9, 10, 11, 12, 13, 14, 15, 20],
    "B0": 3e-2
  }
}
```

### 1.1.2 Example output

```json
{
  "res": {
    "f1": 0.02798768033743801
  }
}
```

## 1.2 Problem description

In this example we compute the uniformity of the magnetic field produced by a distributed winding.

![Fig 1: 2D model from the coil geometry and the design variables](assets/the_problem.jpg)

*Fig 1: 2D model from the coil geometry and the design variables*

This deals with a seemingly simple problem, where the radius of a given number of circular turns should be optimized
that would generate a uniform magnetic field in the prescribed region. The task is to get a region *(control region in Fig. 1)* with a highly uniform magnetic field distribution. This magnetic
field is generated by a prescribed number (N = 10) of massive circular turns of rectangular cross-section *(yellow color
in Fig. 1)*.

The dimensions of the turns and the variation of their positions in the z-direction are fixed. The coil has 10 turns, each with 1 mm x 1.5mm dimension and 3A excitation. The region [0 mm, 5 mm] x [-5 mm , 5 mm] is
called the control region, and we measure the radial and tangential components of the magnetic flux density within this
area.

Each turn in the winding has a particular distance from the z axis. This is the inner radius of the turns (radii) which
can be varied from 5 mm to 50 mm in the r-direction. These distances are the [input of the model](#x-input variable), so ten unknown radii (design variables)
are to be identified.

After the computation we plot the results in a 2D contour plot.

### 1.2.1 Parameters summary

| Parameter                                                                                      | Value                         |
|------------------------------------------------------------------------------------------------|-------------------------------|
| Control region                                                                                 | [0 mm, 5 mm] x [-5 mm , 5 mm] |
| Number of turns                                                                                | 10                            |
| Turn's dimension height (h)                                                                    | 1.5 mm                        |
| Turn's dimension width (w)                                                                     | 1mm                           |
| Turn's excitation                                                                              | 3A                            |
| <span id="radii">Radii (turns distance from z axis) R1,..., Ri,...,R10 of the ten turns</span> | variation range is 5≤Ri≤50 mm |


## 1.3 The $F_1$ function

The goal function for single-objective problem is to design the geometry of coils that minimize discrepancy between the prescribed valued
$\vec{B}_0$ and actual distribution of $\vec{B}$
in the region of interest.

$$
F_1(r) = \sup_{q=1,n_{\rm p}} |\vec{B}(r_q,z_q)-\vec{B}_0(r_q,z_q)|,
$$

where <span id="b0">$B_0(r_q,z_q) = (0, 2\,{\rm mT})$ is the prescribed value of magnetic flux density</span> and
$n_{\rm p}$ is the number of points.

The value of $B$ is calculated in $n_p$ different points of the region of interest.
