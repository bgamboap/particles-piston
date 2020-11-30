#include "gravity_piston.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>

#define PRECISION 9

// compile with
// g++ -Wall -O3 -fopenmp -I/usr/include/eigen3 -DLOG=1 -DGP_DEBUG=0 -DM_DEBUG=0 main.cpp gravity_piston.cpp -o pistonlog

Eigen::Array<double, -1, -1> uniformize(Eigen::Array<double, -1, -1> arr, double Tmax, double dt){
    //std::cout <<"entered av\n"<<std::flush;
    unsigned N = arr.rows();
    unsigned NT = unsigned(Tmax/dt + 0.000001) + 1;
    Eigen::Array<int, -1, -1> first_in_bin(1, NT);
    first_in_bin = -1; // initialize to impossible values
    unsigned LenTest = 7;
    Eigen::Array<double, -1 ,-1> test(LenTest, 2);
    test << 0.0, 1.0, 
         0.01, 1.1, 
         0.02, 1.2, 
         0.54, 2.1, 
         0.55, 2.0, 
         0.61, 1.6,  
         1.02, 0.3;
    //std::cout << "NT:" << NT << "\n";
    //std::cout << test << "\n";
    unsigned bin;
    for(unsigned n = LenTest; n --> 0;){
        //std::cout << "n:"  << n << std::flush;
        bin = test(n,0)/dt;
        //std::cout << " bin:" << bin << "\n" << std::flush;
        first_in_bin(bin) = n;
    }

    //std::cout << "first in bin:\n";
    //std::cout << first_in_bin << "\n";

    for(unsigned n = NT - 1; n --> 0;){
        if(first_in_bin(n) == -1)
            first_in_bin(n) = first_in_bin(n+1);
    }
    //std::cout << "first in bin:\n";
    //std::cout << first_in_bin << "\n";

    Eigen::Array<double, -1, -1> uniform(NT, 2);

    uniform(0,0) = 0;
    uniform(0,1) = test(0,1);
    unsigned n;
    double t, tn, tprev, xn, xprev;
    for(unsigned i = 1; i < NT; i++){
        //std::cout << "i:" << i << std::flush;
        t = i*dt;
        uniform(i,0) = t;
        n = first_in_bin(i);

        //std::cout << " n:" << n << std::flush;
        tn = test(n,0);
        tprev = test(n-1,0);

        xn = test(n,1);
        xprev = test(n-1, 1);

        //std::cout << " i:" << i << std::flush;
        uniform(i,1) = xprev + (xn - xprev)/(tn - tprev)*(t-tprev);
        //std::cout << "\n";
    }

    //std::cout << "uniform:\n";
    //std::cout << uniform << "\n";

    return uniform;
}


int main(int argc, char** argv){

    unsigned NL, NR, NIter, Nmedias;
    double gamma, TL, TR;
    if(argc != 8){
        std::cout << "Wrong number of parameters. Exiting.\n";
        exit(1);
    }

    NL      = atoi(argv[1]);
    NR      = atoi(argv[2]);
    NIter   = atoi(argv[3]);
    gamma   = atof(argv[4]);
    TL      = atof(argv[5]);
    TR      = atof(argv[6]);
    Nmedias = atoi(argv[7]);

    std::cout << "Program parameters:\n";
    std::cout << "NL:      " << NL      << "\n";
    std::cout << "NR:      " << NR      << "\n";
    std::cout << "NIter:   " << NIter   << "\n";
    std::cout << "gamma:   " << gamma   << "\n";
    std::cout << "TL:      " << TL      << "\n";
    std::cout << "TR:      " << TR      << "\n";
    std::cout << "Nmedias: " << Nmedias << "\n";


    //omp_set_num_threads(2);
    double dt = 0.2;
    double Tmax = 1.0;
    //av_iter = Tmax/dt;
    Eigen::Array<double, -1, -1> global_arr(NIter, 5);
#pragma omp parallel
{
#pragma omp for
    for(unsigned n = 0; n < Nmedias; n++){
        Eigen::Array<double, -1, -1> arr, uniform;
        arr = gravity_piston(NL, NR, NIter, gamma, TL, TR, n);
#pragma omp critical
    {
        if(M_DEBUG) std::cout << "id=" << n << " critical\n" << std::flush;
        global_arr += arr;
        uniform = uniformize(arr, Tmax, dt);
    }
    }
#pragma omp barrier
    global_arr /= Nmedias;
}

    // Save the observables to a file
    std::ofstream file2;
    file2.open("observables.dat");
    file2 << std::setprecision(PRECISION) << global_arr;
    file2.close();
    return 0;
}
