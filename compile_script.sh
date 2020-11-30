#!/bin/bash
GP_DEBUG=0
M_DEBUG=0

cd src
 g++ -Wall -O3 -fopenmp -I/usr/include/eigen3 -DM_DEBUG=$M_DEBUG -DGP_DEBUG=$GP_DEBUG -DLOG=1 *.cpp -o pistonlog

 g++ -Wall -O3 -fopenmp -I/usr/include/eigen3 -DM_DEBUG=$M_DEBUG -DGP_DEBUG=$GP_DEBUG -DLOG=0 *.cpp -o piston
 cd ..
