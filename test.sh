#!/bin/sh

readonly PREFIX="set xrange [-5:5]; set yrange [-3:3]; set zrange [-5:5];"
build/msc 3 test.txt | gnuplot -p -e "$PREFIX splot '<cat' using 2:3:4:1 with points palette"
build/test_custom_struct | gnuplot -p -e "$PREFIX splot '<cat' using 2:3:4:1 with points palette"
build/test_1d_flat_vector | gnuplot -p -e "$PREFIX plot '<cat' using 2:1"

