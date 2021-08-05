import agros2d as a2d

# PROBLEM
problem = a2d.problem(clear=True)
problem.coordinate_type = "axisymmetric"
problem.mesh_type = "triangle"

magnetic = a2d.field("magnetic")
magnetic.analysis_type = "steadystate"
magnetic.number_of_refinements = 0
magnetic.polynomial_order = 1
magnetic.solver = "linear"

geometry = a2d.geometry

magnetic.adaptivity_type = "hp-adaptivity"
magnetic.adaptivity_parameters["tolerance"] = 0.0005
magnetic.adaptivity_parameters["steps"] = 100

# MATERIAL DEFINITIONS
magnetic.add_material(
    "J+",
    {
        "magnetic_remanence_angle": 0.0,
        "magnetic_velocity_y": 0.0,
        "magnetic_current_density_external_real": 2000000.0,
        "magnetic_permeability": 1.0,
        "magnetic_conductivity": 0.0,
        "magnetic_remanence": 0.0,
        "magnetic_velocity_angular": 0.0,
        "magnetic_velocity_x": 0.0,
    },
)
magnetic.add_material(
    "air",
    {
        "magnetic_remanence_angle": 0.0,
        "magnetic_velocity_y": 0.0,
        "magnetic_current_density_external_real": 0.0,
        "magnetic_permeability": 1.0,
        "magnetic_conductivity": 0.0,
        "magnetic_remanence": 0.0,
        "magnetic_velocity_angular": 0.0,
        "magnetic_velocity_x": 0.0,
    },
)
magnetic.add_material(
    "control",
    {
        "magnetic_remanence_angle": 0.0,
        "magnetic_velocity_y": 0.0,
        "magnetic_current_density_external_real": 0.0,
        "magnetic_permeability": 1.0,
        "magnetic_conductivity": 0.0,
        "magnetic_remanence": 0.0,
        "magnetic_velocity_angular": 0.0,
        "magnetic_velocity_x": 0.0,
    },
)

# BOUNDARY DEFINITIONS
magnetic.add_boundary("a0", "magnetic_potential", {"magnetic_potential_real": 0.0})

# GEOMETRY
geometry.add_edge(0.07, 0.035, 0.0, 0.035, boundaries={"magnetic": "a0"})
geometry.add_edge(0.0, -0.005, 0.0, -0.035, boundaries={"magnetic": "a0"})
geometry.add_edge(0.0, -0.035, 0.07, -0.035, boundaries={"magnetic": "a0"})
geometry.add_edge(0.0, 0.035, 0.0, 0.005, boundaries={"magnetic": "a0"})
geometry.add_edge(0.0, 0.005, 0.0, -0.005, boundaries={"magnetic": "a0"})
geometry.add_edge(0.07, -0.035, 0.07, 0.035, boundaries={"magnetic": "a0"})
geometry.add_edge(0.0, -0.005, 0.005, -0.005)
geometry.add_edge(0.005, -0.005, 0.005, 0.005)
geometry.add_edge(0.005, 0.005, 0.0, 0.005)
geometry.add_edge(0.0010775719999999989, -0.015, 0.002077572000000003, -0.015)
geometry.add_edge(0.002077572, -0.015, 0.002077572, -0.0135)
geometry.add_edge(0.002077572, -0.0135, 0.0010775720000000012, -0.0135)
geometry.add_edge(0.001077572, -0.0135, 0.001077572, -0.015)
geometry.add_edge(0.012081099, -0.0135, 0.013081099, -0.0135)
geometry.add_edge(0.013081099, -0.0135, 0.013081099, -0.012)
geometry.add_edge(0.013081099, -0.012, 0.012081099, -0.012)
geometry.add_edge(0.012081099, -0.012, 0.012081099, -0.0135)
geometry.add_edge(0.005120719, -0.012, 0.006120719, -0.012)
geometry.add_edge(0.006120719, -0.012, 0.006120719, -0.0105)
geometry.add_edge(0.006120719, -0.0105, 0.005120719, -0.0105)
geometry.add_edge(0.005120719, -0.0105, 0.005120719, -0.012)
geometry.add_edge(0.020172014000000002, -0.0105, 0.021172014000000003, -0.0105)
geometry.add_edge(0.021172014000000003, -0.0105, 0.021172014000000003, -0.009000000000000001)
geometry.add_edge(0.021172014000000003, -0.009000000000000001, 0.020172014000000002, -0.009000000000000001)
geometry.add_edge(0.020172014000000002, -0.009000000000000001, 0.020172014000000002, -0.0105)
geometry.add_edge(0.016271996, -0.009000000000000001, 0.017271996, -0.009000000000000001)
geometry.add_edge(0.017271996, -0.009000000000000001, 0.017271996, -0.0075)
geometry.add_edge(0.017271996, -0.0075, 0.016271996, -0.0075)
geometry.add_edge(0.016271996, -0.0075, 0.016271996, -0.009000000000000001)
geometry.add_edge(0.005526311, -0.0075, 0.006526311, -0.0075)
geometry.add_edge(0.006526311, -0.0075, 0.006526311, -0.006)
geometry.add_edge(0.006526311, -0.006, 0.005526311, -0.006)
geometry.add_edge(0.005526311, -0.006, 0.005526311, -0.0075)
geometry.add_edge(0.02574817800000001, -0.006, 0.026748178, -0.006)
geometry.add_edge(0.026748178, -0.006, 0.026748178, -0.0045000000000000005)
geometry.add_edge(0.026748178, -0.0045000000000000005, 0.025748178, -0.0045000000000000005)
geometry.add_edge(0.025748178, -0.0045000000000000005, 0.025748178, -0.006)
geometry.add_edge(0.009355905, -0.0045000000000000005, 0.010355905, -0.0045000000000000005)
geometry.add_edge(0.010355905, -0.0045000000000000005, 0.010355905, -0.003)
geometry.add_edge(0.010355905, -0.003, 0.009355905, -0.003)
geometry.add_edge(0.009355905, -0.003, 0.009355905, -0.0045000000000000005)
geometry.add_edge(0.045987766, -0.003, 0.046987766, -0.003)
geometry.add_edge(0.046987766, -0.003, 0.046987766, -0.0015)
geometry.add_edge(0.046987766, -0.0015, 0.045987766, -0.0015)
geometry.add_edge(0.045987766, -0.0015, 0.045987766, -0.003)
geometry.add_edge(0.007401189, -0.0015, 0.008401189, -0.0015)
geometry.add_edge(0.008401189, -0.0015, 0.008401189, 0.0)
geometry.add_edge(0.008401189, 0.0, 0.007401189, 0.0)
geometry.add_edge(0.007401189, 0.0, 0.007401189, -0.0015)
geometry.add_edge(0.010443497000000001, 0.0, 0.011443497, 0.0)
geometry.add_edge(0.011443497, 0.0, 0.011443497, 0.0015)
geometry.add_edge(0.011443497, 0.0015, 0.010443497000000001, 0.0015)
geometry.add_edge(0.010443497000000001, 0.0015, 0.010443497000000001, 0.0)
geometry.add_edge(0.038464941, 0.0015, 0.039464941, 0.0015)
geometry.add_edge(0.039464941, 0.0015, 0.039464941, 0.003)
geometry.add_edge(0.039464940999999996, 0.003, 0.038464941, 0.003)
geometry.add_edge(0.038464941, 0.003, 0.038464941, 0.0015)
geometry.add_edge(0.016626341, 0.003, 0.017626341, 0.003)
geometry.add_edge(0.017626341, 0.003, 0.017626341, 0.0045000000000000005)
geometry.add_edge(0.017626341, 0.0045000000000000005, 0.016626341, 0.0045000000000000005)
geometry.add_edge(0.016626341, 0.0045000000000000005, 0.016626341, 0.003)
geometry.add_edge(0.007767936, 0.0045000000000000005, 0.008767936, 0.0045000000000000005)
geometry.add_edge(0.008767936, 0.0045000000000000005, 0.008767936, 0.006)
geometry.add_edge(0.008767936, 0.006, 0.008425435999999996, 0.006)
geometry.add_edge(0.008425435999999996, 0.006, 0.007767936, 0.006)
geometry.add_edge(0.007767936, 0.006, 0.007767936, 0.0045000000000000005)
geometry.add_edge(0.008767936, 0.006, 0.009425435999999999, 0.006)
geometry.add_edge(0.009425435999999999, 0.006, 0.009425435999999999, 0.0075)
geometry.add_edge(0.009425435999999999, 0.0075, 0.008425436, 0.0075)
geometry.add_edge(0.008425436, 0.0075, 0.008425436, 0.006)
geometry.add_edge(0.022641502, 0.0075, 0.023641502, 0.0075)
geometry.add_edge(0.023641502, 0.0075, 0.023641502, 0.009000000000000001)
geometry.add_edge(0.023641502, 0.009000000000000001, 0.022641502, 0.009000000000000001)
geometry.add_edge(0.022641502, 0.009000000000000001, 0.022641502, 0.0075)
geometry.add_edge(0.046928319, 0.009000000000000001, 0.047928319000000004, 0.009000000000000001)
geometry.add_edge(0.047928319000000004, 0.009000000000000001, 0.047928319000000004, 0.0105)
geometry.add_edge(0.047928319000000004, 0.0105, 0.046928319, 0.0105)
geometry.add_edge(0.046928319, 0.0105, 0.046928319, 0.009000000000000001)
geometry.add_edge(0.00569151, 0.0105, 0.00669151, 0.0105)
geometry.add_edge(0.00669151, 0.0105, 0.00669151, 0.012)
geometry.add_edge(0.00669151, 0.012, 0.005915284999999999, 0.012)
geometry.add_edge(0.005915284999999999, 0.012, 0.00569151, 0.012)
geometry.add_edge(0.00569151, 0.012, 0.00569151, 0.0105)
geometry.add_edge(0.004915285, 0.012, 0.005691510000000001, 0.012)
geometry.add_edge(0.005915285, 0.012, 0.005915285, 0.0135)
geometry.add_edge(0.005915285, 0.0135, 0.004915285, 0.0135)
geometry.add_edge(0.004915285, 0.0135, 0.004915285, 0.012)
geometry.add_edge(0.001373434, 0.0135, 0.0023734340000000002, 0.0135)
geometry.add_edge(0.0023734340000000002, 0.0135, 0.0023734340000000002, 0.015)
geometry.add_edge(0.0023734340000000002, 0.015, 0.001373434, 0.015)
geometry.add_edge(0.001373434, 0.015, 0.001373434, 0.0135)

# BLOCK LABELS
geometry.add_label(0.00182757234384367, -0.01425, materials={"magnetic": "J+"})
geometry.add_label(0.012831098619021401, -0.012750000000000001, materials={"magnetic": "J+"})
geometry.add_label(0.00587071874756443, -0.01125, materials={"magnetic": "J+"})
geometry.add_label(0.0209220140666903, -0.00975, materials={"magnetic": "J+"})
geometry.add_label(0.017021996038583002, -0.00825, materials={"magnetic": "J+"})
geometry.add_label(0.00627631142548928, -0.00675, materials={"magnetic": "J+"})
geometry.add_label(0.0264981784448498, -0.00525, materials={"magnetic": "J+"})
geometry.add_label(0.01010590470354692, -0.00375, materials={"magnetic": "J+"})
geometry.add_label(0.0467377661249235, -0.0022500000000000003, materials={"magnetic": "J+"})
geometry.add_label(0.008151189281054209, -0.00075, materials={"magnetic": "J+"})
geometry.add_label(0.0111934974707302, 0.00075, materials={"magnetic": "J+"})
geometry.add_label(0.0392149413574651, 0.0022500000000000003, materials={"magnetic": "J+"})
geometry.add_label(0.0173763411147178, 0.00375, materials={"magnetic": "J+"})
geometry.add_label(0.008517935505111929, 0.00525, materials={"magnetic": "J+"})
geometry.add_label(0.00917543646671979, 0.00675, materials={"magnetic": "J+"})
geometry.add_label(0.023391502259683703, 0.00825, materials={"magnetic": "J+"})
geometry.add_label(0.047678319345423, 0.00975, materials={"magnetic": "J+"})
geometry.add_label(0.0064415104034465905, 0.01125, materials={"magnetic": "J+"})
geometry.add_label(0.0056652852743515605, 0.012750000000000001, materials={"magnetic": "J+"})
geometry.add_label(0.00212343445767236, 0.01425, materials={"magnetic": "J+"})
geometry.add_label(0.03, 0.03, materials={"magnetic": "air"})
geometry.add_label(0.003, 0.0, materials={"magnetic": "control"})

# SOLVE
problem.solve()
a2d.view.zoom_best_fit()
f = open(
    "/home/gadok/work/adze-modeler/applications/Cooking_with_adze_modeler/AgrosvsFemm/snapshots/finest-mesh-solution/agros_solution.csv",
    "w",
)

# POSTPROCESSING AND EXPORTING
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
