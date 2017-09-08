#!/bin/bash

./primDPG $* \
    -file-control    ./files/control \
    -file-geometry   ./files/geometry \
    -file-phys       ./files/physics \
    -file-err        ./files/dump_err \
    -prefix          'pbl1_' \
    -paraview-attr   \
    -vis-level       2 \
    -problem         1 \
    -inner-product   1 \
    -diffusion       1.E-1
