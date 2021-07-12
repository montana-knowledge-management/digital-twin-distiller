
-- PROBLEM

showconsole()
clearconsole()
remove("femm_solution.csv")
newdocument(0)
file_out = openfile("femm_solution.csv", "w")
mi_probdef(0.0,'millimeters','axi',1e-08, 1.0, 30, 0)
smartmesh(0)

-- MATERIAL DEFINITIONS

mi_addmaterial('air', 1.0, 1.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0.0, 0.0)

-- BOUNDARY DEFINITIONS

mi_addboundprop('d0', 30, 30, 30, 0, 0, 0, 0, 0, 0, 0, 0)

-- GEOMETRY

mi_addnode(-1, 0)
mi_addnode(1, 0)
mi_addsegment(-1, 0, 1, 0)

-- BLOCK LABELS

mi_addblocklabel(0, 0)
mi_selectlabel(0, 0)
mi_setblockprop('air', 0, 1.0, '<None>', 0.0, 0, 0)
mi_clearselected()

-- SOLVE

mi_saveas("femm_solver_script.fem")
mi_analyze(1)
mi_loadsolution()

-- POSTPROCESSING AND EXPORTING


-- CLOSING STEPS

closefile(file_out)
mo_close()
mi_close()
quit()
