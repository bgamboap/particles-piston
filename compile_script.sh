#!/bin/bash

# Debug flags. Set to 1 if program has issues
GP_DEBUG=0  # Debug flag in the gravity_piston.cpp file
M_DEBUG=0   # Debug flag in the main.cpp file

cd src
 g++ -Wall -O3 -fopenmp -I/usr/include/eigen3 -DM_DEBUG=$M_DEBUG -DGP_DEBUG=$GP_DEBUG -DLOG=1 *.cpp -o pistonlog

 g++ -Wall -O3 -fopenmp -I/usr/include/eigen3 -DM_DEBUG=$M_DEBUG -DGP_DEBUG=$GP_DEBUG -DLOG=0 *.cpp -o piston
 cd ..
