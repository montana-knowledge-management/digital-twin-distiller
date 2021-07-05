import numpy
import numpy as np

from adze_modeler.boundaries import DirichletBoundaryCondition
from adze_modeler.geometry import Geometry
from adze_modeler.material import Material
from adze_modeler.snapshot import Snapshot
from adze_modeler.platforms.agros2d import Agros2D
from adze_modeler.platforms.femm import Femm
from adze_modeler.metadata import Agros2DMetadata, FemmMetadata
from numpy import linspace, meshgrid, genfromtxt
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy import interpolate
import os

def evaluate(platform, Nx=20, Ny=20):

    snapshot = Snapshot(platform)

    geo = Geometry()
    geo.import_svg("geometry.svg")
    geo.generate_intersections()
    snapshot.add_geometry(geo)

    # boundaries
    b1 = DirichletBoundaryCondition(name="a0", field_type=snapshot.platform.metadata.problem_type,
                                    magnetic_potential=0.0)
    # print(b1)
    snapshot.add_boundary_condition(b1)

    # materials
    exctitation = Material("J+")
    exctitation.Je = 2e6
    air = Material("air")
    core = Material("core")
    core.mu_r = 400
    core.meshsize = 0.1
    snapshot.add_material(exctitation)
    snapshot.add_material(air)
    snapshot.add_material(core)

    snapshot.assign_material(3, 0, name="core")
    snapshot.assign_material(30, 30, name="air")

    # block labeles for coils
    snapshot.assign_material(34.274, -14.250, name='J+')
    snapshot.assign_material(6.625, -12.750, name='J+')
    snapshot.assign_material(17.876, -11.250, name='J+')
    snapshot.assign_material(15.544, -9.750, name='J+')
    snapshot.assign_material(38.641, -8.250, name='J+')
    snapshot.assign_material(35.951, -6.750, name='J+')
    snapshot.assign_material(45.648, -5.250, name='J+')
    snapshot.assign_material(9.412, -3.750, name='J+')
    snapshot.assign_material(24.486, -2.250, name='J+')
    snapshot.assign_material(6.841, -0.750, name='J+')
    snapshot.assign_material(15.339, 0.750, name='J+')
    snapshot.assign_material(28.241, 2.250, name='J+')
    snapshot.assign_material(6.694, 3.750, name='J+')
    snapshot.assign_material(14.448, 5.250, name='J+')
    snapshot.assign_material(34.745, 6.750, name='J+')
    snapshot.assign_material(30.022, 8.250, name='J+')
    snapshot.assign_material(15.420, 9.750, name='J+')
    snapshot.assign_material(32.017, 11.250, name='J+')
    snapshot.assign_material(41.924, 12.750, name='J+')
    snapshot.assign_material(5.792, 14.250, name='J+')

    # assigning boundaries
    snapshot.assign_boundary_condition(0, 0, "a0")
    snapshot.assign_boundary_condition(0, -20, "a0")
    snapshot.assign_boundary_condition(0, 20, "a0")
    snapshot.assign_boundary_condition(35, 35, "a0")
    snapshot.assign_boundary_condition(35, -35, "a0")
    snapshot.assign_boundary_condition(70, 0, "a0")

    # aadding postprocessing steps
    # adding measurement points to the core
    px = linspace(0, 5, Nx)
    py = linspace(-5, 5, Ny)
    xv, yv = meshgrid(px, py, sparse=False, indexing='xy')

    for i in range(Nx):
        for j in range(Ny):
            eval_point = (xv[j, i], yv[j, i])
            snapshot.add_postprocessing('point_value', eval_point, 'Bz')

    snapshot.export()
    snapshot.execute()

def get_core_points(filename, nb_x=200, nb_y=200):
    name, x, y, Bx = genfromtxt(filename, unpack=True, delimiter=',')
    xi = linspace(x.min(), x.max(), nb_x)
    yi = linspace(y.min(), y.max(), nb_y)
    zi = interpolate.griddata((x, y), Bx, (xi[None, :], yi[:, None]), method='cubic')

    return xi, yi, zi


def plot_platform(ax, xi, yi, zi):
    CS = ax.contourf(xi, yi, zi)
    ax.set_title(r"FEMM - B$_z$")
    ax.clabel(CS, inline=1, fontsize=12, fmt='%1.4f', colors='w')


agros_metadata = Agros2DMetadata()
agros_metadata.file_script_name = "agros_solver_script"
agros_metadata.file_metrics_name = "agros_solution.csv"
agros_metadata.problem_type = "magnetic"
agros_metadata.coordinate_type = "axisymmetric"
agros_metadata.analysis_type = "steadystate"
agros_metadata.unit = 1e-3
agros_metadata.nb_refinements = 0
agros_metadata.adaptivity = "hp-adaptivity"
agros_metadata.adaptivity_tol = 0.001


femm_metadata = FemmMetadata()
femm_metadata.problem_type = "magnetic"
femm_metadata.coordinate_type = "axisymmetric"
femm_metadata.file_script_name = "femm_solver_script"
femm_metadata.file_metrics_name = "femm_solution.csv"
femm_metadata.unit = "millimeters"
femm_metadata.smartmesh = False

platform_agros = Agros2D(agros_metadata)
platform_femm = Femm(femm_metadata)

Bgoal = 2e-3
evaluate(platform_agros)
evaluate(platform_femm)

fig, ax = plt.subplots(3, 1, figsize=(12, 20), dpi=220)

xi, yi, agros_zi = get_core_points(agros_metadata.file_metrics_name)
F1_agros = max(numpy.abs(agros_zi-Bgoal).ravel())
plot_platform(ax[0], xi, yi, agros_zi)
ax[0].set_title('Agros2D - Bz')
ax[0].set_xlabel('R [m]')
ax[0].set_ylabel('Z [m]')


xi, yi, femm_zi = get_core_points(femm_metadata.file_metrics_name)
F1_femm = max(numpy.abs(femm_zi-Bgoal).ravel())
plot_platform(ax[1], xi, yi, femm_zi)
ax[1].set_title('FEMM - Bz')
ax[1].set_xlabel('R [mm]')
ax[1].set_ylabel('Z [mm]')


difference = numpy.abs(femm_zi-agros_zi)
CS = ax[2].contourf(xi, yi, difference, cmap='inferno', vmin=0, vmax=difference.max())
ax[2].set_title(r"Difference")
ax[2].clabel(CS, inline=1, fontsize=10, fmt='%1.4e' ,colors='w')
ax[2].set_xlabel('R [mm]')
ax[2].set_ylabel('Z [mm]')

print('Difference')
print("min:", difference.min())
print("max:", difference.max())
print('Fitness')
print("Agros2d F1:", F1_agros)
print("FEMM F1:", F1_femm)

plt.savefig('comparison.png', bbox_inches='tight')

# os.remove(agros_metadata.file_metrics_name)
# os.remove(femm_metadata.file_metrics_name)