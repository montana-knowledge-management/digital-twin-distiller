#!/bin/bash

Xvfb :127 -screen 0 1024x768x16 &
export DISPLAY=:127
wine $FEMM -lua-script=test.lua
cat *.csv
