import json

from numpy import linspace
import requests
from tqdm import tqdm


url = "http://127.0.0.1:5000/process_sim"

post = {
        "simulation": {
            "type": "default",
            "solver": "femm",
            "nb_segments": 100,
            "meshsize": 1.0,
            "adaptivity_tol": 1
            },
        "version": "2021.12"
        }

# Agros
for mi in tqdm(linspace(1, 0.001, 201)):
    post['simulation']['solver'] = 'agros2d'
    post['simulation']['adaptivity_tol'] = mi
    res = requests.post(url, data=json.dumps(post))


# FEMM
for mi in tqdm(linspace(0.2, 0.005, 201)):
    post['simulation']['solver'] = 'femm'
    post['simulation']['meshsize'] = mi
    res = requests.post(url, data=json.dumps(post))
