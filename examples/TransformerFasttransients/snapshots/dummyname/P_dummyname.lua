
-- PROBLEM

showconsole()
clearconsole()
remove("examples/TransformerFasttransients/snapshots/dummyname/S_dummyname.csv")
newdocument(0)
file_out = openfile("examples/TransformerFasttransients/snapshots/dummyname/S_dummyname.csv", "w")
mi_probdef(10000000.0,'meters','planar',1e-12, 1.0, 30, 0)
mi_smartmesh(0)

-- MATERIAL DEFINITIONS

mi_addmaterial('air', 1.0, 1.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
mi_addmaterial('Ii', 1.0, 1.0, 0.0, 0.0625, 0.0, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
mi_addmaterial('Ij', 1.0, 1.0, 0.0, 0.0625, 0.0, 0, 0, 0, 0, 0, 0, 0.0, 0.0)

-- BOUNDARY DEFINITIONS

mi_addboundprop('a0', 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0)

-- GEOMETRY

mi_addnode(0.075, 0.075)
mi_addnode(0.175, 0.075)
mi_addnode(0.175, 0.375)
mi_addnode(0.075, 0.375)
mi_addnode(0.091, 0.343)
mi_addnode(0.095, 0.343)
mi_addnode(0.095, 0.347)
mi_addnode(0.091, 0.347)
mi_addnode(0.091, 0.303)
mi_addnode(0.095, 0.303)
mi_addnode(0.095, 0.307)
mi_addnode(0.091, 0.307)
mi_addsegment(0.075, 0.075, 0.175, 0.075)
mi_selectsegment(0.125, 0.075)
mi_setsegmentprop("a0", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(0.175, 0.375, 0.075, 0.375)
mi_selectsegment(0.125, 0.375)
mi_setsegmentprop("a0", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(0.175, 0.075, 0.175, 0.375)
mi_selectsegment(0.175, 0.225)
mi_setsegmentprop("a0", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(0.075, 0.375, 0.075, 0.075)
mi_selectsegment(0.075, 0.225)
mi_setsegmentprop("a0", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(0.091, 0.343, 0.095, 0.343)
mi_addsegment(0.095, 0.343, 0.095, 0.347)
mi_addsegment(0.095, 0.347, 0.091, 0.347)
mi_addsegment(0.091, 0.347, 0.091, 0.343)
mi_addsegment(0.091, 0.303, 0.095, 0.303)
mi_addsegment(0.095, 0.303, 0.095, 0.307)
mi_addsegment(0.095, 0.307, 0.091, 0.307)
mi_addsegment(0.091, 0.307, 0.091, 0.303)

-- BLOCK LABELS

mi_addblocklabel(0.125, 0.225)
mi_selectlabel(0.125, 0.225)
mi_setblockprop('air', 0, 0, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(0.093, 0.345)
mi_selectlabel(0.093, 0.345)
mi_setblockprop('Ii', 0, 0, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(0.093, 0.305)
mi_selectlabel(0.093, 0.305)
mi_setblockprop('Ij', 0, 0, '<None>', 0.0, 0, 0)
mi_clearselected()

-- SOLVE

mi_saveas("examples/TransformerFasttransients/snapshots/dummyname/P_dummyname.fem")
mi_analyze(1)
mi_loadsolution()

-- POSTPROCESSING AND EXPORTING

mo_selectblock(0.125, 0.225)
Energy = mo_blockintegral(2)
mo_clearblock()
write(file_out, "Energy, ", Energy, "\n")
write(file_out, "nodes, ", mo_numnodes(), "\n")
write(file_out, "elements, ", mo_numelements(), "\n")

-- CLOSING STEPS

closefile(file_out)
mo_close()
mi_close()
quit()
