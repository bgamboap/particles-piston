#!/bin/bash

 g++ -Wall -O3 -fopenmp -I/usr/include/eigen3 -DLOG=1 gravity_piston.cpp -o pistonlog

 g++ -Wall -O3 -fopenmp -I/usr/include/eigen3 -DLOG=0 gravity_piston.cpp -o piston
