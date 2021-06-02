from adze_modeler.dxf_handlers import import_dxf
from importlib_resources import files

geo = import_dxf("2horse.dxf")

if __name__ == "__main__":
    print(geo)
