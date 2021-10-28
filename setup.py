from setuptools import setup

# Metadata goes in setup.cfg. These are here for GitHub's dependency graph.
setup(
    name="Adze-Modeler",
    install_requires=[
        "pre-commit",
        "numpy",
        "svgpathtools",
        "importlib-resources",
        "pygmsh",
        "gmsh",
        "ezdxf",
        "networkx",
        "h5py",
        "scipy",
        "pyvista",
        "vedo",
        "aiofiles",
        "jinja2",
        "fastapi",
        "pydantic",
        "uvicorn",
    ],
    extras_require={"full": ["agrossuite >= 0.01"]},
)
