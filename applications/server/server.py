import uvicorn
from pathlib import Path
from fastapi import FastAPI
from typing import Optional
import sys
import multiprocessing
from numpy import linspace

DIR_APPS = Path(__file__).parent.parent
adze_server = FastAPI()

bldcpath = f'{DIR_APPS.resolve()}'
sys.path.append(bldcpath)
from BLDC.model import BLDCMotor

@adze_server.get('/models/BLDC')
async def model_bldc(rotorangle: Optional[float]=0.0,
                     alpha: Optional[float]=0.0,
                     I0: Optional[float]=0.0,
                     ):
    m = BLDCMotor(rotorangle=rotorangle, alpha=alpha, I0=I0)

    return m(devmode=False, cleanup=True)

def execute_model(m: BLDCMotor):
    return m(cleanup=True).get('Torque', -1996) * 8

@adze_server.get('/models/BLDC/cogging')
async def model_bldc(theta0: Optional[float]=0.0,
                     theta1: Optional[float]=360/24,
                     nb_steps: Optional[int]=11,
                    ):
    theta = linspace(theta0, theta1, nb_steps)
    models = [BLDCMotor(rotorangle=ti) for ti in theta]
    with multiprocessing.Pool(processes=4) as pool:
        T = pool.map(execute_model, models)

    return {'theta': list(theta), 'Torque': T}

if __name__ == "__main__":
    uvicorn.run(adze_server, host="127.0.0.1", port=8000)
