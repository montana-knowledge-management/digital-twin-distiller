import multiprocessing
import numpy as np
import csv
from model import PriusMotor


def execute_model(model_i):
    T = model_i(cleanup=True, devmode=False).get('Torque')
    print(model_i.rotorangle, T)
    return T

def calculate_torque():
    theta = np.linspace(0, 1, 201) * 360
    T = []

    # Generating models
    models = [PriusMotor(rotorangle=theta_i) for theta_i in theta]
    
    # Executing the models
    with multiprocessing.Pool(processes=4) as pool:
        T = pool.map(execute_model, models)

    # Export results
    with open(models[0].dir_data / "torque0.csv", "w") as f:
        for theta_i, Ti in zip(theta, T):
            f.write(f'{theta_i}, {Ti}\n')


def read_results(filename):

    T = []
    theta = []
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            row = [float(ri) for ri in row]
            theta.append(row[0])
            T.append(row[1])

    return theta, T

def curve_fit(x, y, order=3):
    return np.poly1d(np.polyfit(x, y, order))
    
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    m = PriusMotor(rotorangle=0, exportname='pathtest')
    # print(m(cleanup=False, devmode=False))
    theta, T = read_results(m.dir_data/'torque0.csv')
    T1 = curve_fit(theta, T, order=15)

    theta_fine = np.linspace(0, 1, 1001) * 360
    plt.figure(figsize=(6,4))
    plt.plot(theta_fine, T1(theta_fine), 'b-')
    plt.grid()
    plt.xlim(45, 45+180)
    plt.xlabel('Rotor angle [°]')
    plt.ylabel('Torque [Nm]')
    plt.show()

