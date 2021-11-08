
[![tests](https://github.com/robust-design-stack/adze-modeler/actions/workflows/ci.yml/badge.svg)](https://github.com/robust/actions)
[![codecov](https://codecov.io/gh/montana-knowledge-management/digital-twin-distiller/branch/main/graph/badge.svg?token=FPRAPGB6AY)](https://codecov.io/gh/montana-knowledge-management/digital-twin-distiller)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/robust-design-stack/adze-modeler.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/robust-design-stack/adze-modeler/alerts/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/robust-design-stack/adze-modeler.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/robust-design-stack/adze-modeler/context:python)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Code Style](https://badgen.net/badge/Code%20Style/black?labelColor=2e3a44&color=000000)](https://github.com/psf/black)
<img alt="GitHub commit activity" src="https://img.shields.io/github/commit-activity/m/robust-design-stack/Adze-modeler">


## INSTALL

> pip3 install adze-modeler

### For Developers
Digital-Twin-Distiller project uses poetry for dependency management

Installing poetry on your system

The pip package is deprecated, please download teh installation file
> curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py | python -

Run the installation script
> python3 install-poetry.py

Adding poetry to the system path add this line to .bashrc
> export PATH="$HOME/.poetry/bin:$PATH"

If not compiling
>rm -rf ~/.cache/pypoetry

### Run tests and calculating the coverage automatically

> poetry run coverage run

Generating an html from the results
> poetry run coverage html

Install commit-msg git hook. It runs on every local commit to check if the commit message conforms to the convention specified in .gitlint

>pre-commit install --hook-type commit-msg --overwrite
>pre-commit install --hook-type=pre-commit --overwrite
