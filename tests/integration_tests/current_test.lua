showconsole()
clearconsole()
remove("current_data.csv")
newdocument(3)
file_out = openfile("current_data.csv", "w")
ci_probdef("millimeters", "planar", 0, 1e-08, 5, 30)
ci_addnode(-19.0, 5.0)
ci_addnode(19.0, 5.0)
ci_addnode(-19.0, -5.0)
ci_addnode(19.0, -5.0)
ci_addnode(-19.0, 18.0)
ci_addnode(19.0, 18.0)
ci_addnode(-19.0, -18.0)
ci_addnode(19.0, -18.0)
ci_addnode(-42.5, 18.0)
ci_addnode(42.5, 18.0)
ci_addnode(-42.5, -18.0)
ci_addnode(42.5, -18.0)
ci_addsegment(-42.5, 18.0, -42.5, -18.0)
ci_addsegment(-19.0, 18.0, -19.0, 5.0)
ci_addsegment(-19.0, 18.0, -19.0, -5.0)
ci_addsegment(-19.0, -5.0, -19.0, -18.0)
ci_addsegment(42.5, 18.0, 42.5, -18.0)
ci_addsegment(19.0, 18.0, 19.0, 5.0)
ci_addsegment(19.0, 18.0, 19.0, -5.0)
ci_addsegment(19.0, -5.0, 19.0, -18.0)
ci_addsegment(-42.5, 18.0, -19.0, 18.0)
ci_addsegment(-19.0, 18.0, 19.0, 18.0)
ci_addsegment(19.0, 18.0, 42.5, 18.0)
ci_addsegment(-42.5, -18.0, -19.0, -18.0)
ci_addsegment(-19.0, -18.0, 19.0, -18.0)
ci_addsegment(19.0, -18.0, 42.5, -18.0)
ci_addsegment(-19.0, 5.0, 19.0, 5.0)
ci_addsegment(-19.0, -5.0, 19.0, -5.0)
ci_addnode(-28.3, 7.8999999999999995)
ci_addnode(-36.5, 7.8999999999999995)
ci_addnode(-28.3, -7.8999999999999995)
ci_addnode(-36.5, -7.8999999999999995)
ci_addnode(28.3, 7.8999999999999995)
ci_addnode(36.5, 7.8999999999999995)
ci_addnode(28.3, -7.8999999999999995)
ci_addnode(36.5, -7.8999999999999995)
ci_addarc(-28.3, 7.8999999999999995, -36.5, 7.8999999999999995, 180, 1)
ci_addarc(-36.5, 7.8999999999999995, -28.3, 7.8999999999999995, 180, 1)
ci_addarc(-28.3, -7.8999999999999995, -36.5, -7.8999999999999995, 180, 1)
ci_addarc(-36.5, -7.8999999999999995, -28.3, -7.8999999999999995, 180, 1)
ci_addarc(28.3, 7.8999999999999995, 36.5, 7.8999999999999995, 180, 1)
ci_addarc(36.5, 7.8999999999999995, 28.3, 7.8999999999999995, 180, 1)
ci_addarc(28.3, -7.8999999999999995, 36.5, -7.8999999999999995, 180, 1)
ci_addarc(36.5, -7.8999999999999995, 28.3, -7.8999999999999995, 180, 1)
ci_addblocklabel(30.75, 0)
ci_addblocklabel(-30.75, 0)
ci_addblocklabel(0, 11.5)
ci_addblocklabel(0, -11.5)
ci_addblocklabel(32.4, 7.8999999999999995)
ci_addblocklabel(32.4, -7.8999999999999995)
ci_addblocklabel(-32.4, 7.8999999999999995)
ci_addblocklabel(-32.4, -7.8999999999999995)
ci_addblocklabel(0, 0)
ci_addmaterial("Copper", 58000000.0, 58000000.0, 0, 0, 0, 0)
ci_addmaterial("Copper-Manganin", 20833000.0, 20833000.0, 0, 0, 0, 0)
ci_addmaterial("Titanium", 1789000.0, 1789000.0, 0, 0, 0, 0)
ci_selectlabel(30.75, 0)
ci_selectlabel(-30.75, 0)
ci_setblockprop("Copper", 1, 1, 0)
ci_clearselected()
ci_selectlabel(0, 11.5)
ci_selectlabel(0, -11.5)
ci_setblockprop("Copper-Manganin", 1, 1, 0)
ci_clearselected()
ci_selectlabel(32.4, 7.8999999999999995)
ci_selectlabel(32.4, -7.8999999999999995)
ci_selectlabel(-32.4, 7.8999999999999995)
ci_selectlabel(-32.4, -7.8999999999999995)
ci_setblockprop("Titanium", 1, 1, 0)
ci_clearselected()
ci_selectlabel(0, 0)
ci_setblockprop("<No Mesh>", 1, 1, 0)
ci_clearselected()
ci_addboundprop("Jin", 0, 2222222.2222222225, 0, 0, 2)
ci_addboundprop("GND", 0, 0, 0, 0, 0)
ci_selectsegment(-42.5, 0)
ci_setsegmentprop("Jin", 1, 1, 0, 0, "<None>")
ci_clearselected()
ci_selectsegment(42.5, 0)
ci_setsegmentprop("GND", 1, 1, 0, 0, "<None>")
ci_clearselected()
ci_zoomnatural()
ci_zoomout()
hideconsole()
ci_saveas("current_test.fec")
ci_analyze(1)
ci_loadsolution()
co_selectblock(30.75, 0)
co_selectblock(-30.75, 0)
co_selectblock(0, 11.5)
co_selectblock(0, -11.5)
co_selectblock(32.4, 7.8999999999999995)
co_selectblock(32.4, -7.8999999999999995)
co_selectblock(-32.4, 7.8999999999999995)
co_selectblock(-32.4, -7.8999999999999995)
P = co_blockintegral(0)
write(file_out, 'P', ', ', P, "\n") 

closefile(file_out)
quit()