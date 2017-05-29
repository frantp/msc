#!/bin/sh

readonly PREFIX="set datafile separator ','; set xrange [-5:5]; set yrange [-3:3]; set zrange [-5:5];"
build/msc 3 test.csv | gnuplot -p -e "$PREFIX splot '<cat' using 2:3:4:1 with points palette"

