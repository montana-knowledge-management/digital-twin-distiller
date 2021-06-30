import agros2d as a2d

# PROBLEM
problem = a2d.problem(clear=True)
problem.coordinate_type = "axisymmetric"
problem.mesh_type = "triangle"

magnetic = a2d.field("magnetic")
magnetic.analysis_type = "steadystate"
magnetic.number_of_refinements = 2
magnetic.polynomial_order = 2
magnetic.adaptivity_type = "disabled"
magnetic.solver = "linear"

geometry = a2d.geometry


# MATERIAL DEFINITIONS
magnetic.add_material("air", {'magnetic_remanence_angle': 0.0, 'magnetic_velocity_y': 0.0, 'magnetic_current_density_external_real': 0.0, 'magnetic_permeability': 1.0, 'magnetic_conductivity': 0.0, 'magnetic_remanence': 0.0, 'magnetic_velocity_angular': 0.0, 'magnetic_velocity_x': 0.0})

# BOUNDARY DEFINITIONS
magnetic.add_boundary("d0", "magnetic_potential", {'magnetic_potential_real': 30})

# GEOMETRY
geometry.add_edge(-0.001, 0.0, 0.001, 0.0)

# BLOCK LABELS
geometry.add_label(0.0, 0.0, materials = {'magnetic' : 'air'})

# SOLVE
problem.solve()
a2d.view.zoom_best_fit()
f = open("agros_solution.csv", "w")

# POSTPROCESSING AND EXPORTING

# CLOSING STEPS
f.close()
