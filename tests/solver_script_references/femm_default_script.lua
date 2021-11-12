
-- PROBLEM

remove("")
newdocument(0)
file_out = openfile("", "w")
mi_probdef(0.0,'millimeters','axi',1e-08, 1.0, 30, 0)
mi_smartmesh(0)

-- MATERIAL DEFINITIONS

mi_addmaterial('air', 1.0, 1.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)

-- BOUNDARY DEFINITIONS

mi_addboundprop('d0', 30, 30, 30, 0, 0, 0, 0, 0, 0, 0, 0)
mi_addboundprop('cekla', 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0)
mi_addboundprop('retek', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
mi_addboundprop('mogyoro', 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
mi_addboundprop('g', 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0)
mi_addboundprop('f0', 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0)

-- GEOMETRY

mi_addnode(-1, 0)
mi_addnode(1, 0)
mi_addsegment(-1, 0, 1, 0)

-- BLOCK LABELS

mi_addblocklabel(0.0, 0.0)
mi_selectlabel(0.0, 0.0)
mi_setblockprop('air', 0, 0, '<None>', 0.0, 0, 0)
mi_clearselected()

-- SOLVE

mi_saveas("")
mi_analyze(1)
mi_loadsolution()

-- POSTPROCESSING AND EXPORTING


-- CLOSING STEPS

closefile(file_out)
mo_close()
mi_close()
quit()
