# import sys
#
# import adze_modeler.geometry as geo
# import adze_modeler.objects as obj
# import ezdxf
#
#
# def import_dxf(dxf_file):
#     try:
#         doc = ezdxf.readfile(dxf_file)
#     except OSError:
#         print("Not a DXF file or a generic I/O error.")
#         sys.exit(1)
#     except ezdxf.DXFStructureError:
#         print("Invalid or corrupted DXF file.")
#         sys.exit(2)
#
#     # iterate over all entities in modelspace
#     imported_geo = geo.Geometry()
#
#     # id start from the given number
#     id = 0
#
#     msp = doc.modelspace()
#     for e in msp:
#         if e.dxftype() == "LINE":
#             start = geo.Node(e.dxf.start[0], e.dxf.start[1], id)
#             end = geo.Node(e.dxf.end[0], e.dxf.end[1], id + 1)
#             imported_geo.add_line(geo.Line(start, end, id + 2))
#             id += 3
#
#         if e.dxftype() == "ARC":
#             start = geo.Node(e.dxf.start[0], e.dxf.start[1], id)
#             end = geo.Node(e.dxf.end[0], e.dxf.end[1], id + 1)
#             center = geo.Node(e.dxf.center[0], e.dxf.ceter[0], id + 2)
#             imported_geo.add_arc(geo.CircleArc(start, center, end, id + 3))
#             id += 4
#
#         if e.dxftype() == "POLYLINE":
#             print(e.__dict__)
#     return imported_geo
