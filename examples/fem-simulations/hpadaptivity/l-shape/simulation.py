import json
from math import hypot
from uuid import uuid4

import pymongo

from digital_twin_distiller.encapsulator import Encapsulator
from digital_twin_distiller.modelpaths import ModelDir
from digital_twin_distiller.simulationproject import sim
from model import LShape

def post_results(params, r):
    params.pop('type')
    if params['solver'] == 'agros2d':
        params.pop('meshsize')
    else:
        params.pop('adaptivity_tol')

    post = params.copy()

    V = r.pop('V')
    post['V0'] = {'x': V[0][0], 'y':V[0][1], 'V': V[0][2]}
    post['V1'] = {'x': V[1][0], 'y':V[1][1], 'V': V[1][2]}
    post['V2'] = {'x': V[2][0], 'y':V[2][1], 'V': V[2][2]}


    Ex, Ey = r.pop('Ex'), r.pop('Ey')
    post['E0'] = {'x':Ex[0][0], 'y': Ex[0][1], 'Ex': Ex[0][2], 'Ey': Ey[0][2], 'abs': hypot(Ex[0][2], Ey[0][2])}
    post['E1'] = {'x':Ex[1][0], 'y': Ex[1][1], 'Ex': Ex[1][2], 'Ey': Ey[1][2], 'abs': hypot(Ex[1][2], Ey[1][2])}


    post['elements'] = int(r.pop('elements'))
    post['Energy'] = r.pop('Energy')

    conn_str = "mongodb://localhost:27017"
    client = pymongo.MongoClient(conn_str, serverSelectionTimeoutMS=5000, w=0)

    try:
        db = client.lshape
        db.results.insert_one(post.copy())
        post['dbupload'] = True

    except Exception:
        post['dbupload'] = False


    return post

@sim.register('default')
def default_simulation(model, modelparams, simparams, miscparams):
    m = LShape(**simparams)
    res_ = m(cleanup=True, devmode=False, timeout=1000)
    p = post_results(simparams, res_)

    # with open(ModelDir.DATA / f"{uuid4()}.json", "w", encoding="utf-8") as f:
    #     json.dump(p, f)

    return p


if __name__ == "__main__":

    ModelDir.set_base(__file__)

    # set the model for the simulation
    sim.set_model(LShape)

    model = Encapsulator(sim)
    model.build_docs()
    model.run()
