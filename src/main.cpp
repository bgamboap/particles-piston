#include "gravity_piston.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>

#define PRECISION 9

// compile with
// g++ -Wall -O3 -fopenmp -I/usr/include/eigen3 -DLOG=1 -DGP_DEBUG=0 -DM_DEBUG=0 main.cpp gravity_piston.cpp -o pistonlog

int main(int argc, char** argv){

    unsigned NL, NR, NIter, Nmedias;
    double gamma, betat_R, betat_L;
    if(argc != 8){
        std::cout << "Wrong number of parameters. Exiting.\n";
        exit(1);
    }

    NL      = atoi(argv[1]); // number of particles in the left chamber of the piston
    NR      = atoi(argv[2]); // number of particles in the right chamber of the piston
    NIter   = atoi(argv[3]); // number of iterations we want, or number of collisions
    gamma   = atof(argv[4]); // mass of one particle divided by the mass of the wall
    betat_L = atof(argv[5]); // initial temperature of the left chamber
    betat_R = atof(argv[6]); // initial temperature of the right chamber
    Nmedias = atoi(argv[7]); // unique identifier assigned by the paralelization

    std::cout << "Program parameters:\n";
    std::cout << "NL:      " << NL      << "\n";
    std::cout << "NR:      " << NR      << "\n";
    std::cout << "NIter:   " << NIter   << "\n";
    std::cout << "gamma:   " << gamma   << "\n";
    std::cout << "betatil_L:      " << betat_L      << "\n";
    std::cout << "betatil_R:      " << betat_R      << "\n";
    std::cout << "Nmedias: " << Nmedias << "\n";


    unsigned index = 0;
#pragma omp parallel
{
#pragma omp for
    for(unsigned n = 0; n < Nmedias; n++){
        Eigen::Array<double, -1, -1> arr;
        arr = gravity_piston(NL, NR, NIter, gamma, betat_L, betat_R, n);
#pragma omp critical
    {
        if(M_DEBUG) std::cout << "id=" << n << " critical\n" << std::flush;

        // Save the observables to a file
        std::ofstream file2;
        std::string filename;
        filename = std::string("observables") + std::to_string(index) + ".dat";
        file2.open(filename);
        file2 << std::setprecision(PRECISION) << arr;
        file2.close();
        index++;
    }
    }
}
    return 0;
}
