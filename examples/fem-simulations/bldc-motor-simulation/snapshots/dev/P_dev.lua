
-- PROBLEM

remove("/home/oroszt/Projects/digital-twin-distiller/examples/bldc-motor-simulation/snapshots/dev/S_dev.csv")
newdocument(0)
file_out = openfile("/home/oroszt/Projects/digital-twin-distiller/examples/bldc-motor-simulation/snapshots/dev/S_dev.csv", "w")
mi_probdef(0.0,'millimeters','planar',1e-08, 50.0, 30, 0)
mi_smartmesh(0)

-- MATERIAL DEFINITIONS

mi_addmaterial('stator_steel', 1.0, 1.0, 0.0, 0.0, 1.9, 0.635, 0.0, 0.98, 0, 0, 0, 0.0, 1.0)
mi_addbhpoint("stator_steel", 0.0, 0.0)
mi_addbhpoint("stator_steel", 0.05, 15.120714)
mi_addbhpoint("stator_steel", 0.1, 22.718292)
mi_addbhpoint("stator_steel", 0.15, 27.842733)
mi_addbhpoint("stator_steel", 0.2, 31.871434)
mi_addbhpoint("stator_steel", 0.25, 35.365044)
mi_addbhpoint("stator_steel", 0.3, 38.600588)
mi_addbhpoint("stator_steel", 0.35, 41.736202)
mi_addbhpoint("stator_steel", 0.4, 44.873979)
mi_addbhpoint("stator_steel", 0.45, 48.087807)
mi_addbhpoint("stator_steel", 0.5, 51.437236)
mi_addbhpoint("stator_steel", 0.55, 54.975221)
mi_addbhpoint("stator_steel", 0.6, 58.752993)
mi_addbhpoint("stator_steel", 0.65, 62.823644)
mi_addbhpoint("stator_steel", 0.7, 67.245285)
mi_addbhpoint("stator_steel", 0.75, 72.084406)
mi_addbhpoint("stator_steel", 0.8, 77.4201)
mi_addbhpoint("stator_steel", 0.85, 83.350021)
mi_addbhpoint("stator_steel", 0.9, 89.999612)
mi_addbhpoint("stator_steel", 0.95, 97.537353)
mi_addbhpoint("stator_steel", 1.0, 106.201406)
mi_addbhpoint("stator_steel", 1.05, 116.348464)
mi_addbhpoint("stator_steel", 1.1, 128.547329)
mi_addbhpoint("stator_steel", 1.15, 143.765431)
mi_addbhpoint("stator_steel", 1.2, 163.754169)
mi_addbhpoint("stator_steel", 1.25, 191.868158)
mi_addbhpoint("stator_steel", 1.3, 234.833507)
mi_addbhpoint("stator_steel", 1.35, 306.509769)
mi_addbhpoint("stator_steel", 1.4, 435.255202)
mi_addbhpoint("stator_steel", 1.45, 674.911968)
mi_addbhpoint("stator_steel", 1.5, 1108.325569)
mi_addbhpoint("stator_steel", 1.55, 1813.085468)
mi_addbhpoint("stator_steel", 1.6, 2801.217421)
mi_addbhpoint("stator_steel", 1.65, 4053.653117)
mi_addbhpoint("stator_steel", 1.7, 5591.10689)
mi_addbhpoint("stator_steel", 1.75, 7448.318413)
mi_addbhpoint("stator_steel", 1.8, 9708.81567)
mi_addbhpoint("stator_steel", 1.85, 12486.931615)
mi_addbhpoint("stator_steel", 1.9, 16041.483644)
mi_addbhpoint("stator_steel", 1.95, 21249.420624)
mi_addbhpoint("stator_steel", 2.0, 31313.495878)
mi_addbhpoint("stator_steel", 2.05, 53589.446877)
mi_addbhpoint("stator_steel", 2.1, 88477.484601)
mi_addbhpoint("stator_steel", 2.15, 124329.41054)
mi_addbhpoint("stator_steel", 2.2, 159968.5693)
mi_addbhpoint("stator_steel", 2.25, 197751.604272)
mi_addbhpoint("stator_steel", 2.3, 234024.751347)
mi_addmaterial('U+', 1.0, 1.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)
mi_addmaterial('U-', 1.0, 1.0, 0.0, -0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)
mi_addmaterial('V+', 1.0, 1.0, 0.0, -0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)
mi_addmaterial('V-', 1.0, 1.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)
mi_addmaterial('W+', 1.0, 1.0, 0.0, -0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)
mi_addmaterial('W-', 1.0, 1.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)
mi_addmaterial('air', 1.0, 1.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)
mi_addmaterial('airgap', 1.0, 1.0, 0.0, 0.0, 0.0, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)
mi_addmaterial('magnet', 1.11, 1.11, 724000, 0.0, 1.176, 0, 0.0, 0, 0, 0, 0, 0.0, 1.0)
mi_addmaterial('rotor_steel', 1.0, 1.0, 0.0, 0.0, 5.8, 0, 20, 0, 0, 0, 0, 0.0, 1.0)
mi_addbhpoint("rotor_steel", 0.0, 0.0)
mi_addbhpoint("rotor_steel", 0.2503, 238.7325)
mi_addbhpoint("rotor_steel", 0.925, 795.775)
mi_addbhpoint("rotor_steel", 1.25, 1591.55)
mi_addbhpoint("rotor_steel", 1.39, 2387.325)
mi_addbhpoint("rotor_steel", 1.525, 3978.875)
mi_addbhpoint("rotor_steel", 1.71, 7957.75)
mi_addbhpoint("rotor_steel", 1.87, 15915.5)
mi_addbhpoint("rotor_steel", 1.955, 23873.25)
mi_addbhpoint("rotor_steel", 2.02, 39788.75)
mi_addbhpoint("rotor_steel", 2.11, 79577.5)
mi_addbhpoint("rotor_steel", 2.225, 159155.0)
mi_addbhpoint("rotor_steel", 2.43, 318310.0)

-- BOUNDARY DEFINITIONS

mi_addboundprop('a0', 0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0)
mi_addboundprop('PB1', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
mi_addboundprop('PB2', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
mi_addboundprop('PB3', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
mi_addboundprop('PB4', 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0)
mi_addboundprop('APairgap', 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 3.75)

-- GEOMETRY

mi_addnode(-4.362591128962023, 10.532226670628669)
mi_addnode(-9.662756667218515, 23.32795819590999)
mi_addnode(4.362591128962023, 10.532226670628669)
mi_addnode(9.662756667218515, 23.32795819590999)
mi_addnode(-7.928300000000001, 23.97299645663846)
mi_addnode(7.928300000000001, 23.97299645663846)
mi_addnode(-7.9282999999999975, 26.384551523761022)
mi_addnode(7.9282999999999975, 26.384551523761022)
mi_addnode(-10.61946524813124, 25.637657027188208)
mi_addnode(10.61946524813124, 25.637657027188208)
mi_addnode(-10.734270277840768, 25.914820886941595)
mi_addnode(-10.810806964313786, 26.09959679344385)
mi_addnode(-19.134171618254488, 46.19397662556434)
mi_addnode(10.734270277840768, 25.914820886941595)
mi_addnode(10.810806964313786, 26.09959679344385)
mi_addnode(19.134171618254488, 46.19397662556434)
mi_addnode(6.570855, 27.475195)
mi_addnode(6.752029, 28.151343)
mi_addnode(5.807299, 28.78826)
mi_addnode(7.413365, 40.98756)
mi_addnode(8.04707, 27.079645)
mi_addnode(8.228243, 27.755793)
mi_addnode(9.364862, 27.835014)
mi_addnode(14.073617, 39.202951)
mi_addnode(12.995308, 41.509326)
mi_addnode(9.500396, 42.445785)
mi_addnode(-0.764145, 28.239663)
mi_addnode(-0.764145, 28.939663)
mi_addnode(-1.84153, 29.310363)
mi_addnode(-3.4476, 41.509663)
mi_addnode(0.764145, 28.239663)
mi_addnode(0.764145, 28.939663)
mi_addnode(1.84153, 29.310364)
mi_addnode(3.4476, 41.509663)
mi_addnode(1.809099, 43.458363)
mi_addnode(-1.8091, 43.458363)
mi_addnode(-8.04707, 27.079644)
mi_addnode(-8.228243, 27.755793)
mi_addnode(-9.364862, 27.835014)
mi_addnode(-14.073618, 39.202951)
mi_addnode(-6.570855, 27.475195)
mi_addnode(-6.752029, 28.151343)
mi_addnode(-5.807299, 28.78826)
mi_addnode(-7.413366, 40.98756)
mi_addnode(-9.500396, 42.445785)
mi_addnode(-12.995308, 41.509326)
mi_addarc(19.134171618254488, 46.19397662556434, -19.134171618254488, 46.19397662556434, 45.0, 1)
mi_selectarcsegment(0.0, 50.0)
mi_setarcsegmentprop(1, 'a0', 0, 0)
mi_clearselected()
mi_addarc(4.362591128962023, 10.532226670628669, -4.362591128962023, 10.532226670628669, 45.0, 1)
mi_selectarcsegment(0.0, 11.4)
mi_setarcsegmentprop(1, 'a0', 0, 0)
mi_clearselected()
mi_addsegment(-10.810806964313786, 26.09959679344385, -19.134171618254488, 46.19397662556434)
mi_selectsegment(-14.972489291284138, 36.1467867095041)
mi_setsegmentprop("PB1", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(10.810806964313786, 26.09959679344385, 19.134171618254488, 46.19397662556434)
mi_selectsegment(14.972489291284138, 36.1467867095041)
mi_setsegmentprop("PB1", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(10.734270277840768, 25.914820886941595, 10.810806964313786, 26.09959679344385)
mi_selectsegment(10.772538621077278, 26.007208840192725)
mi_setsegmentprop("PB2", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(-10.734270277840768, 25.914820886941595, -10.810806964313786, 26.09959679344385)
mi_selectsegment(-10.772538621077278, 26.007208840192725)
mi_setsegmentprop("PB2", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(9.662756667218515, 23.32795819590999, 10.61946524813124, 25.637657027188208)
mi_selectsegment(10.141110957674877, 24.4828076115491)
mi_setsegmentprop("PB3", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(-9.662756667218515, 23.32795819590999, -10.61946524813124, 25.637657027188208)
mi_selectsegment(-10.141110957674877, 24.4828076115491)
mi_setsegmentprop("PB3", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(-4.362591128962023, 10.532226670628669, -9.662756667218515, 23.32795819590999)
mi_selectsegment(-7.012673898090269, 16.93009243326933)
mi_setsegmentprop("PB4", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addsegment(4.362591128962023, 10.532226670628669, 9.662756667218515, 23.32795819590999)
mi_selectsegment(7.012673898090269, 16.93009243326933)
mi_setsegmentprop("PB4", None, 1, 0, 0, "<None>")
mi_clearselected()
mi_addarc(10.734270277840768, 25.914820886941595, -10.734270277840768, 25.914820886941595, 45.0, 1)
mi_selectarcsegment(0.0, 28.05)
mi_setarcsegmentprop(1, 'APairgap', 0, 0)
mi_clearselected()
mi_addarc(10.61946524813124, 25.637657027188208, -10.61946524813124, 25.637657027188208, 45.0, 1)
mi_selectarcsegment(0.0, 27.75)
mi_setarcsegmentprop(1, 'APairgap', 0, 0)
mi_clearselected()
mi_addsegment(-7.928300000000001, 23.97299645663846, 7.928300000000001, 23.97299645663846)
mi_addsegment(-7.928300000000001, 23.97299645663846, -7.9282999999999975, 26.384551523761022)
mi_addsegment(7.928300000000001, 23.97299645663846, 7.9282999999999975, 26.384551523761022)
mi_addsegment(6.570855, 27.475195, 6.752029, 28.151343)
mi_addsegment(6.752029, 28.151343, 5.807299, 28.78826)
mi_addsegment(5.807299, 28.78826, 7.413365, 40.98756)
mi_addsegment(8.04707, 27.079645, 8.228243, 27.755793)
mi_addsegment(8.228243, 27.755793, 9.364862, 27.835014)
mi_addsegment(9.364862, 27.835014, 14.073617, 39.202951)
mi_addsegment(5.807299, 28.78826, 9.364862, 27.835014)
mi_addsegment(-0.764145, 28.239663, -0.764145, 28.939663)
mi_addsegment(-0.764145, 28.939663, -1.84153, 29.310363)
mi_addsegment(-1.84153, 29.310363, -3.4476, 41.509663)
mi_addsegment(0.764145, 28.239663, 0.764145, 28.939663)
mi_addsegment(0.764145, 28.939663, 1.84153, 29.310364)
mi_addsegment(1.84153, 29.310364, 3.4476, 41.509663)
mi_addsegment(-1.84153, 29.310363, 1.84153, 29.310364)
mi_addsegment(-8.04707, 27.079644, -8.228243, 27.755793)
mi_addsegment(-8.228243, 27.755793, -9.364862, 27.835014)
mi_addsegment(-9.364862, 27.835014, -14.073618, 39.202951)
mi_addsegment(-6.570855, 27.475195, -6.752029, 28.151343)
mi_addsegment(-6.752029, 28.151343, -5.807299, 28.78826)
mi_addsegment(-5.807299, 28.78826, -7.413366, 40.98756)
mi_addsegment(-9.364862, 27.835014, -5.807299, 28.78826)
mi_addarc(-7.928300000000001, 23.97299645663846, -9.662756667218515, 23.32795819590999, 4.2, 1)
mi_addarc(9.662756667218515, 23.32795819590999, 7.928300000000001, 23.97299645663846, 4.2, 1)
mi_addarc(7.9282999999999975, 26.384551523761022, -7.9282999999999975, 26.384551523761022, 33.45, 1)
mi_addarc(10.810806964313786, 26.09959679344385, -10.810806964313786, 26.09959679344385, 45.0, 1)
mi_addarc(12.995308, 41.509326, 9.500396, 42.445785, 4.77, 20)
mi_addarc(9.500396, 42.445785, 7.413365, 40.98756, 95.12, 20)
mi_addarc(14.073617, 39.202951, 12.995308, 41.509326, 95.12, 20)
mi_addarc(1.809099, 43.458363, -1.8091, 43.458363, 4.77, 20)
mi_addarc(-1.8091, 43.458363, -3.4476, 41.509663, 95.12, 20)
mi_addarc(3.4476, 41.509663, 1.809099, 43.458363, 95.12, 20)
mi_addarc(-9.500396, 42.445785, -12.995308, 41.509326, 4.77, 20)
mi_addarc(-12.995308, 41.509326, -14.073618, 39.202951, 95.12, 20)
mi_addarc(-7.413366, 40.98756, -9.500396, 42.445785, 95.12, 20)

-- BLOCK LABELS

mi_addblocklabel(0.0, 49.9)
mi_selectlabel(0.0, 49.9)
mi_setblockprop('stator_steel', 0, 1.2, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(9.940458, 33.995606)
mi_selectlabel(9.940458, 33.995606)
mi_setblockprop('U+', 0, 1.0, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(0.803035, 35.410014)
mi_selectlabel(0.803035, 35.410014)
mi_setblockprop('V-', 0, 1.0, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(-8.389114, 34.411287)
mi_selectlabel(-8.389114, 34.411287)
mi_setblockprop('W+', 0, 1.0, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(7.399549, 27.615494)
mi_selectlabel(7.399549, 27.615494)
mi_setblockprop('air', 0, 1.0, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(-0.0, 28.589663)
mi_selectlabel(-0.0, 28.589663)
mi_setblockprop('air', 0, 1.0, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(-7.39955, 27.615494)
mi_selectlabel(-7.39955, 27.615494)
mi_setblockprop('air', 0, 1.0, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(-8.795528333609257, 24.856254859835506)
mi_selectlabel(-8.795528333609257, 24.856254859835506)
mi_setblockprop('airgap', 0, 0.18, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(0.0, 28.15)
mi_selectlabel(0.0, 28.15)
mi_setblockprop('airgap', 0, 0.18, '<None>', 0.0, 0, 0)
mi_clearselected()
mi_addblocklabel(0.0, 26.4)
mi_selectlabel(0.0, 26.4)
mi_setblockprop('magnet', 0, 0.18, '<None>', 90.0, 0, 0)
mi_clearselected()
mi_addblocklabel(0.0, 18.325)
mi_selectlabel(0.0, 18.325)
mi_setblockprop('rotor_steel', 0, 0.18, '<None>', 0.0, 0, 0)
mi_clearselected()

-- SOLVE

mi_saveas("/home/oroszt/Projects/digital-twin-distiller/examples/bldc-motor-simulation/snapshots/dev/P_dev.fem")
mi_analyze(1)
mi_loadsolution()

-- POSTPROCESSING AND EXPORTING

mo_selectblock(1.1220826297187624e-15, 18.325)
mo_selectblock(1.6165337748745062e-15, 26.4)
Torque = mo_blockintegral(22)
mo_clearblock()
write(file_out, "Torque, ", Torque, "\n")

-- CLOSING STEPS

closefile(file_out)
mo_close()
mi_close()
quit()
