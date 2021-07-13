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
geometry.add_edge(8.8817842e-19, 0.14000000000000004, 8.8817842e-19, 0.004999999903142452, boundaries={'magnetic': 'a0'})
geometry.add_edge(8.8817842e-19, 0.004999999903142452, 8.8817842e-19, 2.842170943040401e-17, boundaries={'magnetic': 'a0'})
geometry.add_edge(0.14, 2.8421709e-17, 0.14, 0.14000000000000004, boundaries={'magnetic': 'a0'})
geometry.add_edge(0.14000000000000004, 0.14000000000000004, 0.0, 0.14000000000000004, boundaries={'magnetic': 'a0'})
geometry.add_edge(0.018485976999999997, 2.8421709e-17, 0.14, 2.8421709e-17, boundaries={'magnetic': 'n0'})
geometry.add_edge(0.005000000000000001, 2.8421709e-17, 0.017485977, 2.8421709e-17, boundaries={'magnetic': 'n0'})
geometry.add_edge(0.017485977, 2.8421709e-17, 0.018485976999999997, 2.8421709e-17, boundaries={'magnetic': 'n0'})
geometry.add_edge(1.9721522630525296e-34, 2.8421709e-17, 0.005000000000000001, 2.8421709e-17, boundaries={'magnetic': 'n0'})
geometry.add_edge(0.005, 2.842170899365265e-17, 0.005, 0.004999999903142452)
geometry.add_edge(0.005, 0.004999999903142452, 8.881784197001253e-19, 0.004999999903142452)
geometry.add_edge(0.018485976999999997, 2.8421708999999994e-17, 0.018485976999999997, 0.0015)
geometry.add_edge(0.018485976999999997, 0.0015, 0.017485977, 0.0015)
geometry.add_edge(0.017485977, 0.0015, 0.017485977, 2.842170943040401e-17)
geometry.add_edge(0.014401638000000001, 0.0015, 0.015401638, 0.0015)
geometry.add_edge(0.015401638, 0.0015, 0.015401638, 0.003)
geometry.add_edge(0.015401638, 0.003, 0.014427641000000003, 0.003)
geometry.add_edge(0.014427641000000003, 0.003, 0.014401638000000001, 0.003)
geometry.add_edge(0.014401638000000001, 0.003, 0.014401638000000001, 0.0015)
geometry.add_edge(0.013427641, 0.003, 0.014401638000000003, 0.003)
geometry.add_edge(0.014427641, 0.003, 0.014427641, 0.0045000000000000005)
geometry.add_edge(0.014427641, 0.0045000000000000005, 0.013427641, 0.0045000000000000005)
geometry.add_edge(0.013427641, 0.0045000000000000005, 0.013427641, 0.003)
geometry.add_edge(0.019451646, 0.0045000000000000005, 0.020451646, 0.0045000000000000005)
geometry.add_edge(0.020451646, 0.0045000000000000005, 0.020451646, 0.006)
geometry.add_edge(0.020451646, 0.006, 0.019451646, 0.006)
geometry.add_edge(0.019451646, 0.006, 0.019451646, 0.0045000000000000005)
geometry.add_edge(0.002668015, 0.006, 0.003668015, 0.006)
geometry.add_edge(0.003668015, 0.006, 0.003668015, 0.0075)
geometry.add_edge(0.0036680150000000023, 0.0075, 0.002668015, 0.0075)
geometry.add_edge(0.002668015, 0.0075, 0.002668015, 0.006)
geometry.add_edge(0.004936238999999999, 0.0075, 0.005936238999999999, 0.0075)
geometry.add_edge(0.005936238999999999, 0.0075, 0.005936238999999999, 0.009000000000000001)
geometry.add_edge(0.005936238999999999, 0.009000000000000001, 0.004936238999999999, 0.009000000000000001)
geometry.add_edge(0.004936238999999999, 0.009000000000000001, 0.004936238999999999, 0.0075)
geometry.add_edge(0.019909513, 0.009000000000000001, 0.020909513, 0.009000000000000001)
geometry.add_edge(0.020909513, 0.009000000000000001, 0.020909513, 0.0105)
geometry.add_edge(0.020909513, 0.0105, 0.019909513, 0.0105)
geometry.add_edge(0.019909513, 0.0105, 0.019909513, 0.009000000000000001)
geometry.add_edge(0.007290303, 0.0105, 0.008290303, 0.0105)
geometry.add_edge(0.008290303, 0.0105, 0.008290303, 0.012)
geometry.add_edge(0.008290303, 0.012, 0.008162454999999999, 0.012)
geometry.add_edge(0.008162454999999999, 0.012, 0.007290303, 0.012)
geometry.add_edge(0.007290303, 0.012, 0.007290303, 0.0105)
geometry.add_edge(0.007162455, 0.012, 0.007290303, 0.012)
geometry.add_edge(0.008162454999999999, 0.012, 0.008162454999999999, 0.0135)
geometry.add_edge(0.008162454999999999, 0.0135, 0.007162455, 0.0135)
geometry.add_edge(0.007162455, 0.0135, 0.007162455, 0.012)
geometry.add_edge(0.008766665, 0.0135, 0.009766665, 0.0135)
geometry.add_edge(0.009766665, 0.0135, 0.009766665, 0.015)
geometry.add_edge(0.009766665, 0.015, 0.008766665, 0.015)
geometry.add_edge(0.008766665, 0.015, 0.008766665, 0.0135)

# BLOCK LABELS
geometry.add_label(0.018235977040213, 0.00075, materials = {'magnetic' : 'J+'})
geometry.add_label(0.015151638466662001, 0.0022500000000000003, materials = {'magnetic' : 'J+'})
geometry.add_label(0.014177640833697, 0.00375, materials = {'magnetic' : 'J+'})
geometry.add_label(0.020201645879223, 0.00525, materials = {'magnetic' : 'J+'})
geometry.add_label(0.0034180148538810002, 0.00675, materials = {'magnetic' : 'J+'})
geometry.add_label(0.005686238693079, 0.00825, materials = {'magnetic' : 'J+'})
geometry.add_label(0.020659512969778, 0.00975, materials = {'magnetic' : 'J+'})
geometry.add_label(0.008040302720339, 0.01125, materials = {'magnetic' : 'J+'})
geometry.add_label(0.007912454911862, 0.012750000000000001, materials = {'magnetic' : 'J+'})
geometry.add_label(0.009516665045609, 0.01425, materials = {'magnetic' : 'J+'})
geometry.add_label(0.03, 0.03, materials = {'magnetic' : 'air'})
geometry.add_label(0.003, 0.001, materials = {'magnetic' : 'control'})

# SOLVE
problem.solve()
a2d.view.zoom_best_fit()
f = open(r"/home/gadokrisztian/work/adze-modeler/applications/Cooking_with_adze_modeler/ArtapOptimizationSymmetricNewF1/snapshots/72be41b3-cea5-4dd0-9271-94724366b55e/agros_solution.csv", "w")

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
