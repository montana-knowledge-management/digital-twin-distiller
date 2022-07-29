import sys

from digital_twin_distiller.cli import optimize_cli

if __name__ == "__main__":
    args = sys.argv.copy()
    args.pop(0)
    optimize_cli(args)
