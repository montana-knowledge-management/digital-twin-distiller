import argparse
import string
import subprocess
from os import chdir, getcwd
from pathlib import Path
from shutil import copy

from digital_twin_distiller.modelpaths import ModelDir


def new(name, location):
    """
    Creates a project template in the given location under the given name: ~/location/name
    :parameter name: creates a new project with the given name
    :parameter location: creates a project under the given location
    """
    location = Path(location)
    if location.suffix:
        location = location.parent

    SRC = Path(__file__).parent.resolve() / "resources"
    SRC_CODE = SRC / "model_template"
    SRC_DOC = SRC / "doc_template"
    DST = Path(location).resolve() / name
    ModelDir.set_base(DST)

    # Creating the directory tree
    for dir_i in ModelDir.get_dirs():
        # print(dir_i)
        dir_i.mkdir(exist_ok=True, parents=True)

    # copy template files
    for file_i in SRC_CODE.iterdir():
        copy(file_i, DST / file_i.name)

    # copy the docs template
    for file_i in SRC_DOC.rglob("*"):
        if not file_i.is_dir():
            folder = file_i.relative_to(SRC_DOC).parent
            fname = file_i.name

            dst = DST / "docs" / folder
            dst.mkdir(exist_ok=True, parents=True)

            copy(file_i, dst / fname)

    #  default json-s
    with open(ModelDir.DEFAULTS / "model.json", "w") as f:
        f.write("{}")

    with open(ModelDir.DEFAULTS / "simulation.json", "w") as f:
        print("{", '  "default": {}', "}", sep="\n", file=f)

    with open(ModelDir.DEFAULTS / "misc.json", "w") as f:
        f.write("{}")

    # replace the model name in the files
    for file_i in DST.rglob("*"):
        if not file_i.is_dir() and file_i.suffix in {".py", ".md", ".yml"}:
            with open(file_i, encoding="utf-8") as f:
                template = string.Template(f.read())

            with open(file_i, "w", encoding="utf-8") as f:
                f.write(template.substitute(name=name))

    # build the documentation
    cwd = getcwd()
    chdir(DST / "docs")
    subprocess.run(["mkdocs", "build", "-q"])
    chdir(cwd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="ADZE-modeler", description="Create a new Model")
    parser.add_argument("task", help="Select a task: [new, ]")
    parser.add_argument("name", help="The name of the model", default="MODEL")
    parser.add_argument("location", help="The location of the model", default="APPLICATIONS")
    args = parser.parse_args()

    if args.task == "new":
        new(args.name, args.location)
