#!/bin/bash

./primDPG $* \
    -file-control    ./files/control \
    -file-geometry   ./files/geometry \
    -file-phys       ./files/physics \
    -file-err        ./files/dump_err \
    -prefix          'pbl2_' \
    -paraview-attr   \
    -vis-level       2 \
    -problem         2 \
    -inner-product   2 \
    -diffusion       1.E-1
