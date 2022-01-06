
[![tests](https://github.com/robust-design-stack/adze-modeler/actions/workflows/ci.yml/badge.svg)](https://github.com/robust/actions)
[![codecov](https://codecov.io/gh/montana-knowledge-management/digital-twin-distiller/branch/main/graph/badge.svg?token=FPRAPGB6AY)](https://codecov.io/gh/montana-knowledge-management/digital-twin-distiller)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/robust-design-stack/adze-modeler.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/robust-design-stack/adze-modeler/alerts/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/robust-design-stack/adze-modeler.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/robust-design-stack/adze-modeler/context:python)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Code Style](https://badgen.net/badge/Code%20Style/black?labelColor=2e3a44&color=000000)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
<img alt="GitHub commit activity" src="https://img.shields.io/github/commit-activity/m/robust-design-stack/Adze-modeler">

# Digital-Twin-Distiller

Python project for creating a long-lasting, encapsulated version of your numerical simulation or your machine-learning-based project.

## INSTALL

For manage the packages in a Linux operating system use the `apt` package management tool.

With `apt list` you can see all the available packages.

Before trying to install a package you should update the tool:

```shell
sudo apt-get update #optional
```

1. If Python 3 has already deployed (check by asking for the version: `python3 --version`) install **pip3**:

```shell
sudo apt-get -y install python3-pip
```

*Verify pip3 installation: `pip3 --version`*

2. Install the project

```shell
pip3 install digital-twin-distiller
```

### For Developers

Digital-Twin-Distiller project uses **poetry** for dependency management.

#### Installing poetry on your system

1. Install python-is-python3 package

```shell
sudo apt-get install -y python-is-python3
```

*Once the installation is done, python command will use Python 3.x binary.*

*Verification: `python --version`.*

2. Install venv package

For running *install-poetry.py* script we need the *venv* package to be installed.

On Debian/Ubuntu systems, you need to install the *python3-venv* package for creating the virtual environment.

```shell
sudo apt install -y python3.8-venv
```

3. The pip package is deprecated, please download the installation file
```shell
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py |  python -
```

4. Adding poetry to the system path add this line to *.bashrc* file (`nano .bashrc`)

```shell
export PATH="$HOME/.poetry/bin:$PATH"
```

5. Reloading .bashrc configuration

```shell
. ~/.bashrc
```

6. Verify

```shell
poetry --version
```

7. If not compiling

```shell
rm -rf ~/.cache/pypoetry
```


### Run tests and calculating the coverage automatically

> Clone the project and step into the main folder (digital-twin-distiller).

1. Add the coverage tool

```shell
poetry add -D coverage[toml]
```

2. Calculating coverage

```shell
poetry run coverage run
```

3. Generating a html from the results

```shell
poetry run coverage report
```

4. Install pre-commit package manager to enable running hooks using poetry

```shell
poetry run pre-commit install
```

5. Install commit-msg git hook. It runs on every local commit to check if the commit message conforms to the convention specified in .gitlint

```shell
pre-commit install --hook-type commit-msg --overwrite

pre-commit install --hook-type=pre-commit --overwrite
```
