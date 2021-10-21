#!/bin/bash

export DISPLAY=:127
Xvfb :127 -screen 0 1024x768x16 &

python3 simulation.py
