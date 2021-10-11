"""
how to create a new model:
    python -m adze_modeler.newmodel 
"""
import argparse
from os import getcwd
from pathlib import Path
from adze_modeler.modelpaths import ModelDir
import string

parser = argparse.ArgumentParser(description='Create a new Model')
parser.add_argument("-n", "--name", help="The name of the model", default="MODEL")
parser.add_argument("-d", "--dir", help="The base directory of the model", default='applications')
args = parser.parse_args()


ModelDir.set_base(Path(args.dir).resolve() / args.name)

print("Name:", args.name)
print("Base directory:", ModelDir.BASE)

for dir_i in ModelDir.get_dirs():
    # print(ni)
    dir_i.mkdir(exist_ok=True, parents=True)

TEMPLATES = Path(__file__).parent.resolve() /"resources/model_template"
with open(TEMPLATES/'model.py', 'r') as fsrc, open(ModelDir.BASE/'model.py', 'w') as fdst:
    template = string.Template(fsrc.read())
    fdst.write(template.substitute(name=args.name))

with open(TEMPLATES/'simulation.py', 'r') as fsrc, open(ModelDir.BASE/'simulation.py', 'w') as fdst:
    template = string.Template(fsrc.read())
    fdst.write(template.substitute(name=args.name))

with open(ModelDir.DEFAULTS / 'model.json', 'w') as f:
    f.write('{}')

with open(ModelDir.DEFAULTS / 'simulation.json', 'w') as f:
    print('{', '  "default": {}', '}', sep='\n', file=f)

with open(ModelDir.DEFAULTS / 'misc.json', 'w') as f:
    f.write('{}')
