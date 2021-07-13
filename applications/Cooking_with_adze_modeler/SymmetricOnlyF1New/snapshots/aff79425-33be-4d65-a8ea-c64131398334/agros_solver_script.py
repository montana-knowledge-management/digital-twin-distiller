import agros2d as a2d

# PROBLEM
problem = a2d.problem(clear=True)
problem.coordinate_type = "axisymmetric"
problem.mesh_type = "triangle"

magnetic = a2d.field("magnetic")
magnetic.analysis_type = "steadystate"
magnetic.number_of_refinements = 0
magnetic.polynomial_order = 2
magnetic.solver = "linear"

geometry = a2d.geometry

magnetic.adaptivity_type = "hp-adaptivity"
magnetic.adaptivity_parameters["tolerance"] = 1
magnetic.adaptivity_parameters["steps"] = 10

# MATERIAL DEFINITIONS
magnetic.add_material("J+", {'magnetic_remanence_angle': 0.0, 'magnetic_velocity_y': 0.0, 'magnetic_current_density_external_real': 2000000.0, 'magnetic_permeability': 1.0, 'magnetic_conductivity': 0.0, 'magnetic_remanence': 0.0, 'magnetic_velocity_angular': 0.0, 'magnetic_velocity_x': 0.0})
magnetic.add_material("air", {'magnetic_remanence_angle': 0.0, 'magnetic_velocity_y': 0.0, 'magnetic_current_density_external_real': 0.0, 'magnetic_permeability': 1.0, 'magnetic_conductivity': 0.0, 'magnetic_remanence': 0.0, 'magnetic_velocity_angular': 0.0, 'magnetic_velocity_x': 0.0})
magnetic.add_material("control", {'magnetic_remanence_angle': 0.0, 'magnetic_velocity_y': 0.0, 'magnetic_current_density_external_real': 0.0, 'magnetic_permeability': 1.0, 'magnetic_conductivity': 0.0, 'magnetic_remanence': 0.0, 'magnetic_velocity_angular': 0.0, 'magnetic_velocity_x': 0.0})

# BOUNDARY DEFINITIONS
magnetic.add_boundary("a0", "magnetic_potential", {'magnetic_potential_real': 0.0})
magnetic.add_boundary("n0", "magnetic_surface_current", {'magnetic_surface_current_real': 0.0})

# GEOMETRY
geometry.add_edge(8.8817842e-19, 0.004999999903142452, 8.8817842e-19, 2.842170943040401e-17, boundaries={'magnetic': 'a0'})
geometry.add_edge(0.14, 2.8421709e-17, 0.14, 0.14000000000000004, boundaries={'magnetic': 'a0'})
geometry.add_edge(8.8817842e-19, 0.14000000000000004, 8.8817842e-19, 0.004999999903142452, boundaries={'magnetic': 'a0'})
geometry.add_edge(0.14000000000000004, 0.14000000000000004, 0.0, 0.14000000000000004, boundaries={'magnetic': 'a0'})
geometry.add_edge(0.010500312, 2.8421709e-17, 0.14, 2.8421709e-17, boundaries={'magnetic': 'n0'})
geometry.add_edge(0.009500311999999999, 2.8421709e-17, 0.010500312, 2.8421709e-17, boundaries={'magnetic': 'n0'})
geometry.add_edge(1.9721522630525296e-34, 2.8421709e-17, 0.005000000000000001, 2.8421709e-17, boundaries={'magnetic': 'n0'})
geometry.add_edge(0.005000000000000001, 2.8421709e-17, 0.009500311999999997, 2.8421709e-17, boundaries={'magnetic': 'n0'})
geometry.add_edge(0.005, 2.842170899365265e-17, 0.005, 0.004999999903142452)
geometry.add_edge(0.005, 0.004999999903142452, 8.881784197001253e-19, 0.004999999903142452)
geometry.add_edge(0.010500312, 2.8421708999999994e-17, 0.010500312, 0.0015)
geometry.add_edge(0.010500312, 0.0015, 0.009500311999999999, 0.0015)
geometry.add_edge(0.009500311999999999, 0.0015, 0.009500311999999999, 2.842170943040401e-17)
geometry.add_edge(0.013526132000000001, 0.0015, 0.014526132, 0.0015)
geometry.add_edge(0.014526132, 0.0015, 0.014526132, 0.003)
geometry.add_edge(0.014526132, 0.003, 0.013526132000000001, 0.003)
geometry.add_edge(0.013526132000000001, 0.003, 0.013526132000000001, 0.0015)
geometry.add_edge(0.018207136, 0.003, 0.019207136, 0.003)
geometry.add_edge(0.019207136, 0.003, 0.019207136, 0.0045000000000000005)
geometry.add_edge(0.019207136000000003, 0.0045000000000000005, 0.018207136, 0.0045000000000000005)
geometry.add_edge(0.018207136, 0.0045000000000000005, 0.018207136, 0.003)
geometry.add_edge(0.015914209, 0.0045000000000000005, 0.016914209, 0.0045000000000000005)
geometry.add_edge(0.016914209, 0.0045000000000000005, 0.016914209, 0.006)
geometry.add_edge(0.016914209, 0.006, 0.015914209, 0.006)
geometry.add_edge(0.015914209, 0.006, 0.015914209, 0.0045000000000000005)
geometry.add_edge(0.010902948, 0.006, 0.011902948, 0.006)
geometry.add_edge(0.011902948, 0.006, 0.011902948, 0.0075)
geometry.add_edge(0.011902948, 0.0075, 0.010902948, 0.0075)
geometry.add_edge(0.010902948, 0.0075, 0.010902948, 0.006)
geometry.add_edge(0.009579065999999999, 0.0075, 0.010579066, 0.0075)
geometry.add_edge(0.010579066, 0.0075, 0.010579066, 0.009000000000000001)
geometry.add_edge(0.010579066, 0.009000000000000001, 0.010468775999999999, 0.009000000000000001)
geometry.add_edge(0.010468775999999999, 0.009000000000000001, 0.009579065999999999, 0.009000000000000001)
geometry.add_edge(0.009579065999999999, 0.009000000000000001, 0.009579065999999999, 0.0075)
geometry.add_edge(0.010579066000000002, 0.009000000000000001, 0.011468776, 0.009000000000000001)
geometry.add_edge(0.011468776, 0.009000000000000001, 0.011468776, 0.0105)
geometry.add_edge(0.011468776, 0.0105, 0.01066541299999999, 0.0105)
geometry.add_edge(0.01066541299999999, 0.0105, 0.010468776, 0.0105)
geometry.add_edge(0.010468776, 0.0105, 0.010468776, 0.009000000000000001)
geometry.add_edge(0.01146877599999999, 0.0105, 0.011665413, 0.0105)
geometry.add_edge(0.011665413, 0.0105, 0.011665413, 0.012)
geometry.add_edge(0.01166541299999999, 0.012, 0.010665412999999999, 0.012)
geometry.add_edge(0.010665412999999999, 0.012, 0.010665412999999999, 0.0105)
geometry.add_edge(0.0032020620000000003, 0.012, 0.0042020619999999995, 0.012)
geometry.add_edge(0.0042020619999999995, 0.012, 0.0042020619999999995, 0.0135)
geometry.add_edge(0.0042020619999999995, 0.0135, 0.0032020620000000003, 0.0135)
geometry.add_edge(0.0032020620000000003, 0.0135, 0.0032020620000000003, 0.012)
geometry.add_edge(0.008846768000000001, 0.0135, 0.009846768, 0.0135)
geometry.add_edge(0.009846768, 0.0135, 0.009846768, 0.015)
geometry.add_edge(0.009846768, 0.015, 0.008846768000000001, 0.015)
geometry.add_edge(0.008846768000000001, 0.015, 0.008846768000000001, 0.0135)

# BLOCK LABELS
geometry.add_label(0.010250311904679, 0.00075, materials = {'magnetic' : 'J+'})
geometry.add_label(0.014276131943163, 0.0022500000000000003, materials = {'magnetic' : 'J+'})
geometry.add_label(0.018957136165879997, 0.00375, materials = {'magnetic' : 'J+'})
geometry.add_label(0.016664208562454, 0.00525, materials = {'magnetic' : 'J+'})
geometry.add_label(0.011652948295462, 0.00675, materials = {'magnetic' : 'J+'})
geometry.add_label(0.010329065700672, 0.00825, materials = {'magnetic' : 'J+'})
geometry.add_label(0.011218775679075, 0.00975, materials = {'magnetic' : 'J+'})
geometry.add_label(0.011415412825745, 0.01125, materials = {'magnetic' : 'J+'})
geometry.add_label(0.003952062184873, 0.012750000000000001, materials = {'magnetic' : 'J+'})
geometry.add_label(0.009596768194978, 0.01425, materials = {'magnetic' : 'J+'})
geometry.add_label(0.03, 0.03, materials = {'magnetic' : 'air'})
geometry.add_label(0.003, 0.001, materials = {'magnetic' : 'control'})

# SOLVE
problem.solve()
a2d.view.zoom_best_fit()
f = open(r"/home/gadokrisztian/work/adze-modeler/applications/Cooking_with_adze_modeler/ArtapOptimizationSymmetricNewF1/snapshots/aff79425-33be-4d65-a8ea-c64131398334/agros_solution.csv", "w")

# POSTPROCESSING AND EXPORTING
point = magnetic.local_values(1e-06, 1e-06)["Brz"]
f.write("{}, 1e-06, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(1e-06, 1e-06)["Brr"]
f.write("{}, 1e-06, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(1e-06, 0.00125075)["Brz"]
f.write("{}, 1e-06, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(1e-06, 0.00125075)["Brr"]
f.write("{}, 1e-06, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(1e-06, 0.0025004999999999997)["Brz"]
f.write("{}, 1e-06, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(1e-06, 0.0025004999999999997)["Brr"]
f.write("{}, 1e-06, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(1e-06, 0.00375025)["Brz"]
f.write("{}, 1e-06, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(1e-06, 0.00375025)["Brr"]
f.write("{}, 1e-06, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(1e-06, 0.005)["Brz"]
f.write("{}, 1e-06, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(1e-06, 0.005)["Brr"]
f.write("{}, 1e-06, 0.005, {}\n".format("Br", point))

point = magnetic.local_values(0.0005564444444444444, 1e-06)["Brz"]
f.write("{}, 0.0005564444444444444, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(0.0005564444444444444, 1e-06)["Brr"]
f.write("{}, 0.0005564444444444444, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(0.0005564444444444444, 0.00125075)["Brz"]
f.write("{}, 0.0005564444444444444, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(0.0005564444444444444, 0.00125075)["Brr"]
f.write("{}, 0.0005564444444444444, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(0.0005564444444444444, 0.0025004999999999997)["Brz"]
f.write("{}, 0.0005564444444444444, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(0.0005564444444444444, 0.0025004999999999997)["Brr"]
f.write("{}, 0.0005564444444444444, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(0.0005564444444444444, 0.00375025)["Brz"]
f.write("{}, 0.0005564444444444444, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(0.0005564444444444444, 0.00375025)["Brr"]
f.write("{}, 0.0005564444444444444, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(0.0005564444444444444, 0.005)["Brz"]
f.write("{}, 0.0005564444444444444, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(0.0005564444444444444, 0.005)["Brr"]
f.write("{}, 0.0005564444444444444, 0.005, {}\n".format("Br", point))

point = magnetic.local_values(0.0011118888888888888, 1e-06)["Brz"]
f.write("{}, 0.0011118888888888888, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(0.0011118888888888888, 1e-06)["Brr"]
f.write("{}, 0.0011118888888888888, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(0.0011118888888888888, 0.00125075)["Brz"]
f.write("{}, 0.0011118888888888888, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(0.0011118888888888888, 0.00125075)["Brr"]
f.write("{}, 0.0011118888888888888, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(0.0011118888888888888, 0.0025004999999999997)["Brz"]
f.write("{}, 0.0011118888888888888, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(0.0011118888888888888, 0.0025004999999999997)["Brr"]
f.write("{}, 0.0011118888888888888, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(0.0011118888888888888, 0.00375025)["Brz"]
f.write("{}, 0.0011118888888888888, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(0.0011118888888888888, 0.00375025)["Brr"]
f.write("{}, 0.0011118888888888888, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(0.0011118888888888888, 0.005)["Brz"]
f.write("{}, 0.0011118888888888888, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(0.0011118888888888888, 0.005)["Brr"]
f.write("{}, 0.0011118888888888888, 0.005, {}\n".format("Br", point))

point = magnetic.local_values(0.0016673333333333332, 1e-06)["Brz"]
f.write("{}, 0.0016673333333333332, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(0.0016673333333333332, 1e-06)["Brr"]
f.write("{}, 0.0016673333333333332, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(0.0016673333333333332, 0.00125075)["Brz"]
f.write("{}, 0.0016673333333333332, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(0.0016673333333333332, 0.00125075)["Brr"]
f.write("{}, 0.0016673333333333332, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(0.0016673333333333332, 0.0025004999999999997)["Brz"]
f.write("{}, 0.0016673333333333332, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(0.0016673333333333332, 0.0025004999999999997)["Brr"]
f.write("{}, 0.0016673333333333332, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(0.0016673333333333332, 0.00375025)["Brz"]
f.write("{}, 0.0016673333333333332, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(0.0016673333333333332, 0.00375025)["Brr"]
f.write("{}, 0.0016673333333333332, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(0.0016673333333333332, 0.005)["Brz"]
f.write("{}, 0.0016673333333333332, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(0.0016673333333333332, 0.005)["Brr"]
f.write("{}, 0.0016673333333333332, 0.005, {}\n".format("Br", point))

point = magnetic.local_values(0.0022227777777777775, 1e-06)["Brz"]
f.write("{}, 0.0022227777777777775, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(0.0022227777777777775, 1e-06)["Brr"]
f.write("{}, 0.0022227777777777775, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(0.0022227777777777775, 0.00125075)["Brz"]
f.write("{}, 0.0022227777777777775, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(0.0022227777777777775, 0.00125075)["Brr"]
f.write("{}, 0.0022227777777777775, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(0.0022227777777777775, 0.0025004999999999997)["Brz"]
f.write("{}, 0.0022227777777777775, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(0.0022227777777777775, 0.0025004999999999997)["Brr"]
f.write("{}, 0.0022227777777777775, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(0.0022227777777777775, 0.00375025)["Brz"]
f.write("{}, 0.0022227777777777775, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(0.0022227777777777775, 0.00375025)["Brr"]
f.write("{}, 0.0022227777777777775, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(0.0022227777777777775, 0.005)["Brz"]
f.write("{}, 0.0022227777777777775, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(0.0022227777777777775, 0.005)["Brr"]
f.write("{}, 0.0022227777777777775, 0.005, {}\n".format("Br", point))

point = magnetic.local_values(0.002778222222222222, 1e-06)["Brz"]
f.write("{}, 0.002778222222222222, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(0.002778222222222222, 1e-06)["Brr"]
f.write("{}, 0.002778222222222222, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(0.002778222222222222, 0.00125075)["Brz"]
f.write("{}, 0.002778222222222222, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(0.002778222222222222, 0.00125075)["Brr"]
f.write("{}, 0.002778222222222222, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(0.002778222222222222, 0.0025004999999999997)["Brz"]
f.write("{}, 0.002778222222222222, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(0.002778222222222222, 0.0025004999999999997)["Brr"]
f.write("{}, 0.002778222222222222, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(0.002778222222222222, 0.00375025)["Brz"]
f.write("{}, 0.002778222222222222, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(0.002778222222222222, 0.00375025)["Brr"]
f.write("{}, 0.002778222222222222, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(0.002778222222222222, 0.005)["Brz"]
f.write("{}, 0.002778222222222222, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(0.002778222222222222, 0.005)["Brr"]
f.write("{}, 0.002778222222222222, 0.005, {}\n".format("Br", point))

point = magnetic.local_values(0.003333666666666666, 1e-06)["Brz"]
f.write("{}, 0.003333666666666666, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(0.003333666666666666, 1e-06)["Brr"]
f.write("{}, 0.003333666666666666, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(0.003333666666666666, 0.00125075)["Brz"]
f.write("{}, 0.003333666666666666, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(0.003333666666666666, 0.00125075)["Brr"]
f.write("{}, 0.003333666666666666, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(0.003333666666666666, 0.0025004999999999997)["Brz"]
f.write("{}, 0.003333666666666666, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(0.003333666666666666, 0.0025004999999999997)["Brr"]
f.write("{}, 0.003333666666666666, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(0.003333666666666666, 0.00375025)["Brz"]
f.write("{}, 0.003333666666666666, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(0.003333666666666666, 0.00375025)["Brr"]
f.write("{}, 0.003333666666666666, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(0.003333666666666666, 0.005)["Brz"]
f.write("{}, 0.003333666666666666, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(0.003333666666666666, 0.005)["Brr"]
f.write("{}, 0.003333666666666666, 0.005, {}\n".format("Br", point))

point = magnetic.local_values(0.003889111111111111, 1e-06)["Brz"]
f.write("{}, 0.003889111111111111, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(0.003889111111111111, 1e-06)["Brr"]
f.write("{}, 0.003889111111111111, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(0.003889111111111111, 0.00125075)["Brz"]
f.write("{}, 0.003889111111111111, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(0.003889111111111111, 0.00125075)["Brr"]
f.write("{}, 0.003889111111111111, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(0.003889111111111111, 0.0025004999999999997)["Brz"]
f.write("{}, 0.003889111111111111, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(0.003889111111111111, 0.0025004999999999997)["Brr"]
f.write("{}, 0.003889111111111111, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(0.003889111111111111, 0.00375025)["Brz"]
f.write("{}, 0.003889111111111111, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(0.003889111111111111, 0.00375025)["Brr"]
f.write("{}, 0.003889111111111111, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(0.003889111111111111, 0.005)["Brz"]
f.write("{}, 0.003889111111111111, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(0.003889111111111111, 0.005)["Brr"]
f.write("{}, 0.003889111111111111, 0.005, {}\n".format("Br", point))

point = magnetic.local_values(0.004444555555555556, 1e-06)["Brz"]
f.write("{}, 0.004444555555555556, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(0.004444555555555556, 1e-06)["Brr"]
f.write("{}, 0.004444555555555556, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(0.004444555555555556, 0.00125075)["Brz"]
f.write("{}, 0.004444555555555556, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(0.004444555555555556, 0.00125075)["Brr"]
f.write("{}, 0.004444555555555556, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(0.004444555555555556, 0.0025004999999999997)["Brz"]
f.write("{}, 0.004444555555555556, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(0.004444555555555556, 0.0025004999999999997)["Brr"]
f.write("{}, 0.004444555555555556, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(0.004444555555555556, 0.00375025)["Brz"]
f.write("{}, 0.004444555555555556, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(0.004444555555555556, 0.00375025)["Brr"]
f.write("{}, 0.004444555555555556, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(0.004444555555555556, 0.005)["Brz"]
f.write("{}, 0.004444555555555556, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(0.004444555555555556, 0.005)["Brr"]
f.write("{}, 0.004444555555555556, 0.005, {}\n".format("Br", point))

point = magnetic.local_values(0.005, 1e-06)["Brz"]
f.write("{}, 0.005, 1e-06, {}\n".format("Bz", point))

point = magnetic.local_values(0.005, 1e-06)["Brr"]
f.write("{}, 0.005, 1e-06, {}\n".format("Br", point))

point = magnetic.local_values(0.005, 0.00125075)["Brz"]
f.write("{}, 0.005, 0.00125075, {}\n".format("Bz", point))

point = magnetic.local_values(0.005, 0.00125075)["Brr"]
f.write("{}, 0.005, 0.00125075, {}\n".format("Br", point))

point = magnetic.local_values(0.005, 0.0025004999999999997)["Brz"]
f.write("{}, 0.005, 0.0025004999999999997, {}\n".format("Bz", point))

point = magnetic.local_values(0.005, 0.0025004999999999997)["Brr"]
f.write("{}, 0.005, 0.0025004999999999997, {}\n".format("Br", point))

point = magnetic.local_values(0.005, 0.00375025)["Brz"]
f.write("{}, 0.005, 0.00375025, {}\n".format("Bz", point))

point = magnetic.local_values(0.005, 0.00375025)["Brr"]
f.write("{}, 0.005, 0.00375025, {}\n".format("Br", point))

point = magnetic.local_values(0.005, 0.005)["Brz"]
f.write("{}, 0.005, 0.005, {}\n".format("Bz", point))

point = magnetic.local_values(0.005, 0.005)["Brr"]
f.write("{}, 0.005, 0.005, {}\n".format("Br", point))

info = magnetic.solution_mesh_info()
f.write("{}, {}\n".format("dofs", info["dofs"]))
f.write("{}, {}\n".format("nodes", info["nodes"]))
f.write("{}, {}\n".format("elements", info["elements"]))

# CLOSING STEPS
f.close()
