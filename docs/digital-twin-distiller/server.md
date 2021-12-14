# Server class

## What Server class is good for?

The Server class creates the API around the defined project. The server currently supports two project types:

* MachineLearningProject
* SimulationProject The API is implemented through FastAPI.

## Endpoints

Currently the Server class creates two different endpoints for the different projects:

* `/process_ml` for MachineLearningProject
* `/process_sim` for SimulationProject

**IMPORTANT:** The Server class currently cannot handle multiple projects at once. Hence, one API would only answer
to **_only one_** of the above mentioned endpoints.
