from model import *
import h5py
import multiprocessing
from numpy import linspace
import matplotlib.pyplot as plt

alpha = linspace(0, 180, 21)
models = [PriusMotor(I0=250, alpha=ai) for ai in alpha]
with multiprocessing.Pool(processes=4) as pool, h5py.File(DIR_DATA / 'datastore.hdf5', "a") as f:
    T = pool.map(execute_model, models)
    try:
        data = f.create_dataset('locked_rotor', (len(T), 2), dtype='f')
    except ValueError:
        f.pop('locked_rotor')
        data = f.create_dataset('locked_rotor', (len(T), 2), dtype='f')
    finally:
        data[...] = list(zip(alpha, T))

    plt.plot(alpha, T, 'b-')
    plt.show()