import argparse
import json
import string
import subprocess
from importlib import metadata
from importlib.metadata import PackageNotFoundError
from os import chdir, getcwd
from pathlib import Path
from shutil import copy
from sys import version_info

from digital_twin_distiller.modelpaths import ModelDir

NAME_OF_THE_PROGRAM = "digital-twin-distiller"

COMMAND_NEW = "new"
COMMAND_NEW_DESC = "Create a new Model"
DEFAULT_MODEL = {
    "x0": 1.0,
    "mw": 5,
}
DEFAULT_SIMULATION = {
    "default": {"t0": 0.0, "t1": 5.3, "nstep": 101},
}
DEFAULT_MISC = {"processes": 4, "cleanup": True}


def optimize_cli(argv=None):
    """
    Create Command line interface and define argument
    """

    parser = argparse.ArgumentParser(
        prog=NAME_OF_THE_PROGRAM,
        formatter_class=argparse.RawTextHelpFormatter,
        prefix_chars="-",
        description=_get_all_metadata(),
        epilog=f"Run {NAME_OF_THE_PROGRAM} COMMAND --help for more information on a command",
    )

    # optional arguments
    parser.add_argument(
        "-v", "--version", action="version", version=_get_version_text(), help="display version information"
    )
    parser.add_argument(
        "-q",
        "--quiet",
        dest="verbose",
        action="store_false",
        default=True,
        help="suppress output",
    )

    subparser = parser.add_subparsers(dest="command")

    #  register new command
    _register_subparser_new(subparser)

    args = parser.parse_args(argv)

    if args.command == COMMAND_NEW:
        new(args.name, args.location)


def _register_subparser_new(subparser):
    parser_new = subparser.add_parser(COMMAND_NEW, help=COMMAND_NEW_DESC, description=COMMAND_NEW_DESC)
    parser_new.add_argument("name", help="The name of the model", default="MODEL")
    parser_new.add_argument("location", help="The location of the model", default="APPLICATIONS")


def _get_version_text():
    __version__ = _get_metadata("Version")
    return "\n".join(
        [
            f"{NAME_OF_THE_PROGRAM} {__version__} \n"
            f"Python {version_info.major}.{version_info.minor}.{version_info.micro}"
        ]
    )


def _get_all_metadata():
    __description__ = _get_metadata("Summary")
    __license__ = _get_metadata("License")
    __author__ = _get_metadata("Author")
    __author_email__ = _get_metadata("Author-email")

    return "\n".join(
        [
            f"Welcome to Digital Twin Distiller!\n",
            f"Description: {__description__} \n ",
            f"Licence: {__license__} \n ",
            f"Authors: {__author__} <{__author_email__}>\n ",
        ]
    )


def _get_metadata(param: str):
    try:
        __mt__ = metadata.metadata(NAME_OF_THE_PROGRAM).get(param)
    except PackageNotFoundError:
        print(f"[tool.poetry] name attribute ({NAME_OF_THE_PROGRAM}) not found in pyproject.toml. It may be changed.")
        __mt__ = "unknown"
        exit(1)
    return __mt__


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
        json.dump(DEFAULT_MODEL, f, indent=2)

    with open(ModelDir.DEFAULTS / "simulation.json", "w") as f:
        json.dump(DEFAULT_SIMULATION, f, indent=2)

    with open(ModelDir.DEFAULTS / "misc.json", "w") as f:
        json.dump(DEFAULT_MISC, f, indent=2)

    # replace the model name in the files
    for file_i in DST.rglob("*"):
        if not file_i.is_dir() and file_i.suffix in {".py", ".md", ".yml"}:
            with open(file_i, encoding="utf-8") as f:
                template = string.Template(f.read())

            with open(file_i, "w", encoding="utf-8") as f:
                f.write(template.substitute(name=name).replace("SimulationModel", name))

    # build the documentation
    cwd = getcwd()
    chdir(DST / "docs")
    subprocess.run(["mkdocs", "build", "-q"])
    chdir(cwd)
