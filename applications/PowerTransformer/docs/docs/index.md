# Short circuit impedance calculation for large Power Transformers

## References


> Orosz, T.; Pánek, D.; Karban, P. FEM Based Preliminary Design Optimization in
> Case of Large Power Transformers. Appl. Sci. 2020, 10, 1361.
> https://doi.org/10.3390/app10041361


> Andersen, O. (1973). Transformer Leakage Flux Program Based on the Finite
> Element Method. IEEE Transactions on Power Apparatus and Systems, PAS-92(2),
> 682–689. doi:10.1109/tpas.1973.293773


> Orosz, T., Karban, P., Pánek, D., & Doležel, I. (2020). FEM-based transformer
> design optimization technique with evolutionary algorithms and geometric
> programming. International Journal of Applied Electromagnetics and Mechanics,
> 1–9. doi:10.3233/jae-209504

.: <https://github.com/tamasorosz/utopya/blob/main/src/optimization_functions_003.py>

## Model inputs

![](images/drawing.svg)

| Name | Value | Unit  | Description                              |
|------|-------|-------|------------------------------------------|
| w1   | 152   | mm    | width of the window                      |
| h1   | 1129  | mm    | height of the window                     |
| r1   | 184   | mm    | inner radius of the window               |
| z1   | 0     | mm    | height of the window                     |
| w2   | 42    | mm    | width of the low voltage coil            |
| h2   | 979   | mm    | height of the low voltage coil           |
| r2   | 198   | mm    | inner radius of the low voltage coil     |
| z2   | 70    | mm    | height of the low voltage coil           |
| w3   | 41    | mm    | width of the high voltage coil           |
| h3   | 979   | mm    | height of the high voltage coil          |
| r3   | 268   | mm    | inner radius of the high voltage coil    |
| z3   | 70    | mm    | height of the high voltage coil          |
| js   | 3.00  | A/mm2 | current density in the low voltage coil  |
| jp   | -3.02 | A/mm2 | current density in the high voltage coil |

## Endpoints

* `/` endpoint, which contains the project documentation, by default, this page
  contains the documentation of the server module. This can be easily replaced
  by a static html page. ADZE-models supports MkDocs (https://www.mkdocs.org/)
  documentation, however, other documentation generator based outputs can be
  inserted to this endpoint.

* `/process` endpoint, where the application waits for an input json as a
  single input.

* `/apidocs` endpoint, which describes the roles and the basic functions of the
  integrated api endpoint as an OpenAPI documentation. Moreover, this html site
  contains a small, integrated **test interface**.

* `/ping` endpoint, to test the functionality and the accessibility of the
  realized server from the client side.

### ADZE-Server stands on the shoulders of giants:

    FastAPI for the web parts.
    Pydantic for the data parts and the JsonSchma validation.
    OpenAPI for the API documentation and testing the basic functionality.
