showconsole()
clearconsole()
remove("heat_data.csv")
newdocument(2)
file_out = openfile("heat_data.csv", "w")
hi_probdef("meters", "planar", 1e-08, 1, 30, "", 0)
hi_addnode(-0.4, -0.2)
hi_addnode(0.4, -0.2)
hi_addnode(0.4, 0.2)
hi_addnode(-0.4, 0.2)
hi_addnode(-0.1, -0.2)
hi_addnode(-0.1, 0.0)
hi_addnode(0.1, 0.0)
hi_addnode(0.1, -0.2)
hi_addsegment(-0.4, -0.2, -0.1, -0.2)
hi_addsegment(-0.1, -0.2, -0.1, 0.0)
hi_addsegment(-0.1, 0.0, 0.1, 0.0)
hi_addsegment(0.1, 0.0, 0.1, -0.2)
hi_addsegment(0.1, -0.2, 0.4, -0.2)
hi_addsegment(0.4, -0.2, 0.4, 0.2)
hi_addsegment(0.4, 0.2, -0.4, 0.2)
hi_addsegment(-0.4, 0.2, -0.4, -0.2)
hi_addmaterial("Material", 1380, 1380, 0, 0)
hi_addblocklabel(0.0, 0.12)
hi_selectlabel(0.0, 0.12)
hi_setblockprop("Material", 1, 1, 0)
hi_addboundprop("ins", 2, 0, 0, 296.15, 0, 0)
hi_addboundprop("conv", 2, 0, 0, 324.15, 50, 0)
hi_addboundprop("fixtemp", 0, 574.15, 0, 0, 0, 0)
hi_selectsegment(-0.4, 0)
hi_selectsegment(0.4, 0)
hi_setsegmentprop("ins", 1, 1, 0, 0, "<None>")
hi_clearselected()
hi_selectsegment(0, 0.2)
hi_setsegmentprop("conv", 1, 1, 0, 0, "<None>")
hi_clearselected()
hi_selectsegment(0.15000000000000002, -0.2)
hi_selectsegment(-0.15000000000000002, -0.2)
hi_selectsegment(-0.1, -0.1)
hi_selectsegment(0.1, -0.1)
hi_selectsegment(0, 0.0)
hi_setsegmentprop("fixtemp", 1, 1, 0, 0, "<None>")
hi_clearselected()
hi_zoomnatural()
hi_zoomout()
hideconsole()
hi_saveas("heatflow_test.feh")
hi_analyze(1)
hi_loadsolution()
ho_selectblock(0, 0.1)
Fx, Fy = ho_blockintegral(3)
write(file_out, 'Fx', ', ', Fx, "\n")

write(file_out, 'Fy', ', ', Fy, "\n")

closefile(file_out)
ho_close()
hi_close()
quit()
