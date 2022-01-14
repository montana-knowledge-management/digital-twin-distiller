import argparse
from importlib import metadata
from importlib.metadata import PackageNotFoundError
from sys import version_info

# from digital_twin_distiller import new
# from digital_twin_distiller.__main__ import task_new

DTD = "Digital-Twin-Distiller"


def optimize_cli():
    """
    Create Command line interface and define argument
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        prefix_chars="-",
        description=_get_all_metadata(),
    )
    parser.add_argument(
        "-v", "--version", action="version", version=_get_version_text(), help="display version information"
    )

    # TODO: add [new] as subcommand ?
    # subparsers = parser.add_subparsers(dest=task_new().prog)
    # subparsers.add_parser(parser.prog)

    parser.parse_args()


def _get_version_text():
    try:
        __version__ = metadata.version(DTD)
    except PackageNotFoundError:
        __version__ = "unknown"
    return "\n".join(
        [f"{DTD} {__version__} \n" f"Python {version_info.major}.{version_info.minor}.{version_info.micro}"]
    )


def _get_metadata(param: str):
    try:
        __mt__ = metadata.metadata(DTD).get(param)
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
            f"Authors: {__author__} <{__author_email__}>\n " f"--------------------------------------------------",
        ]
    )


if __name__ == "__main__":
    optimize_cli()
