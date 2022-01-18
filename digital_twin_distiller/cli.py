import argparse
from importlib import metadata
from importlib.metadata import PackageNotFoundError
from sys import version_info

from digital_twin_distiller.__main__ import new

NAME_OF_THE_PROGRAM = "digital-twin-distiller"

COMMAND_NEW = "new"
COMMAND_NEW_DESC = "Create a new Model"


def optimize_cli():
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

    # positional arguments
    _register_subparser_new(subparser)
    parser.add_argument(COMMAND_NEW, help=COMMAND_NEW_DESC)

    args = parser.parse_args()

    if args.command == COMMAND_NEW:
        new(args.name, args.location)


def _register_subparser_new(subparser):
    parser = subparser.add_parser(COMMAND_NEW)
    parser.description = COMMAND_NEW_DESC
    parser.prog = NAME_OF_THE_PROGRAM

    parser.add_argument("task", help=f"Select a task: [{COMMAND_NEW}, ]")
    parser.add_argument("name", help="The name of the model", default="MODEL")
    parser.add_argument("location", help="The location of the model", default="APPLICATIONS")
    return parser


def _get_version_text():
    try:
        __version__ = metadata.version(NAME_OF_THE_PROGRAM)
    except PackageNotFoundError:
        __version__ = "unknown"
    return "\n".join(
        [
            f"{NAME_OF_THE_PROGRAM} {__version__} \n"
            f"Python {version_info.major}.{version_info.minor}.{version_info.micro}"
        ]
    )


def _get_metadata(param: str):
    try:
        __mt__ = metadata.metadata(NAME_OF_THE_PROGRAM).get(param)
    except PackageNotFoundError:
        __mt__ = "unknown"
    return __mt__


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
