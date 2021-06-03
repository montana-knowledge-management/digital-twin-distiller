from adze_modeler.dxf_handlers import import_dxf

geo = import_dxf("2horse.dxf")

if __name__ == "__main__":
    print(geo)
