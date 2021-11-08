
-- PROBLEM

remove("/home/oroszt/RobustDesignStack/adze-modeler/examples/bldc-motor-simulation/snapshots/dev/S_dev.csv")
newdocument(0)
file_out = openfile("/home/oroszt/RobustDesignStack/adze-modeler/examples/bldc-motor-simulation/snapshots/dev/S_dev.csv", "w")
mi_probdef(0.0,'millimeters','planar',1e-08, 1000, 30, 0)

-- MATERIAL DEFINITIONS

mi_addmaterial('air', 1.0, 1.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)

-- BOUNDARY DEFINITIONS

mi_addboundprop('a0', 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0)

-- GEOMETRY


-- BLOCK LABELS


-- SOLVE

mi_saveas("/home/oroszt/RobustDesignStack/adze-modeler/examples/bldc-motor-simulation/snapshots/dev/P_dev.fem")
mi_analyze(1)
mi_loadsolution()

-- POSTPROCESSING AND EXPORTING

mo_selectblock(0, 0)
Energy = mo_blockintegral(2)
mo_clearblock()
write(file_out, "Energy, ", Energy, "\n")

-- CLOSING STEPS

closefile(file_out)
mo_close()
mi_close()
quit()
