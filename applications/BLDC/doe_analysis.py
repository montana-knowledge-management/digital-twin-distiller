from model import *
from numpy import linspace
from multiprocessing import Pool
from adze_modeler.doe import fullfact
import operator as op
import string
from time import perf_counter
from adze_modeler.utils import csv_write

DIR_SAVE = DIR_DATA / "doe"

class DOEBLDCMotor(BLDCMotor):

    def __init__(self, rotorangle, exportname=None, **kwargs):
        super(DOEBLDCMotor, self).__init__(rotorangle=rotorangle, exportname=exportname)
        self.airgap += kwargs.get('dairgap', 0.0)
        self.r3 += kwargs.get('dmagnet_h', 0.0)
        self.mw += kwargs.get('dmagnet_w', 0.0)
        self.dHc= kwargs.get('dHc', 0.0)
        self.dmur = kwargs.get('dmur', 0.0)


        # upodating the geometry variables
        self.r4 = self.r3+(self.airgap-self.void) / 2 # Rotor + airgap slice
        self.s1 = self.r3 + self.airgap # Stator Inner Radius
        self.s2 = self.s1 + 21.75       # Stator Outer Radius

    def define_materials(self):
        super(DOEBLDCMotor, self).define_materials()
        self.snapshot.materials['magnet'].mu_r += self.dmur
        self.snapshot.materials['magnet'].coercivity += self.dHc

def doe_full_factorial():
    dXnames = ('dairgap', 'dmagnet_h', 'dmagnet_w', 'dHc', 'dmur')
    dXvalues = (0.05, 0.05, 0.05, 5000, 0.05)
    designs = list(fullfact([3]*4))
    with Pool(processes=12) as pool:
        for i, design_i in enumerate(designs):
            fname = f'D-{i:03}.csv'
            design_i = [di-1 for di in design_i]
            disturbances = map(op.mul, design_i, dXvalues)
            dX = {name_i:di for name_i, di in zip(dXnames, disturbances)}
            print(fname, f'{i+1}/{len(designs)} {(i+1)*100/len(designs):.1f} %')

            theta = linspace(0, 360/24/2, 61)
            models = [DOEBLDCMotor(ti, **dX) for ti in theta]
            t0 = perf_counter()
            T = pool.map(execute_model, models)
            t1 = perf_counter()
            
            csv_write(DIR_SAVE / fname, ['rotorangle', 'Torque'], theta, T)
            print(f'\t Calculation time: {t1-t0:.2f} s')


if __name__ == "__main__":
    doe_full_factorial()