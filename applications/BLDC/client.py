import requests
import json

"""
curl -X 'POST' \
  'http://127.0.0.1:5000/process' \
  -H 'accept: application/json' \
  -H 'Content-Type: application/json' \
  -d '{
  "simulation": {
    "type": "basic"
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
  },
  "version": "0.7"
}'
"""

d = {
    "simulation": {
        "type": "tol1",
        "nsteps":3
    },
    "tolerances": {
        "type": "ff",
        "parameters": {
            "r1":0.04,
            "r2": 0.05,
        },
        "variables": ["Torque", "dummy"]
    },
}
print(json.dumps(d))
r = requests.post('http://127.0.0.1:5000/process', data=json.dumps(d))
print(r.json())
