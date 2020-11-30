#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include "gravity_piston.hpp"

#define inf 999999.9
#define KB 1
#define g 10.0

// Verificar edge cases desta implementação
//
// Compilar com 
// g++ -O3 -Wall -fopenmp -I/usr/include/eigen3 -DLOG=0 gravity_piston.cpp -o piston

Eigen::Array<double, -1, -1> gravity_piston(unsigned NL, unsigned NR, unsigned NIter, double gamma, double TL, double TR, unsigned identifier){

    // Get the command line parameters
    Eigen::Array<double, -1, -1> xL(NL,1), xR(NR,1), vL(NL,1), vR(NR,1); //positions and velocities
    Eigen::Array<double, -1, -1> col_timesL(NL,1), col_timesR(NR,1), col_timesPLp(NL,1), col_timesPLm(NL,1), col_timesPRp(NR,1), col_timesPRm(NR,1); // collision times
    double xM, vM; // Instantaneous position and velocity of the piston
    Eigen::Array<double, -1, -1> observables(NIter, 5); // tracking some variables over time, such as xM
    double vM_temp;


    // list of times, positions, velocities and collision times for every particle
#if LOG > 0
    Eigen::Array<double, -1, -1> state_t;
    state_t = Eigen::Array<double, -1, -1>::Zero(NIter, 5*(NL + NR) + 3); // State vector for all iterations
#endif


    // initializations for the velocities and positions
    xM = 0.5;
    vM = 0.;
    xL = (Eigen::Array<double, -1, -1>::Random(NL, 1) + 1.0)/2.0*xM*0.99;
    xR = (Eigen::Array<double, -1, -1>::Random(NR, 1) + 1.0)/2.0*(1-xM)*0.99 + xM;
    vL = Eigen::Array<double, -1, -1>::Random(NL, 1);
    vR = Eigen::Array<double, -1, -1>::Random(NR, 1);

    std::mt19937 rnd;
    std::normal_distribution<double> normal_dist{0., 1. };
    std::random_device r;
    std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
    rnd.seed(seed2);

    for(unsigned i = 0; i < NL; i++){
      	vL(i) = std::sqrt(2 * KB * TL)*normal_dist(rnd);
    } 

    for(unsigned i = 0; i < NR; i++){
      	vR(i) = std::sqrt(2 * KB * TR)*normal_dist(rnd);
    } 

    vM_temp = -0.1;

    double t = 0; // time variable

    int flag = -1;       // 0: left wall col, 1: left piston col, 2: right wall col, 3: right pison col
    unsigned pos = -1;    // position of the particle that will collide

    for(unsigned iter = 0; iter < NIter; iter++){
        if(GP_DEBUG) std::cout << "ITERATION: " << iter << "\n" << std::flush;

        observables(iter,0) = t;
        observables(iter,1) = xM;
        observables(iter,2) = (vL * vL).mean();
        observables(iter,3) = (vR * vR).mean();
        observables(iter,4) = gamma * (vM-vM_temp);

        vM_temp = vM;

        // calculate all the colision times
        col_timesL  = -(xL - 0.0)/vL;       // time until collision with left wall
        col_timesR  =  (1.0 - xR)/vR;       // time until collision with right wall
        //col_timesPL =  (xM - xL)/(vL - vM); // time until collision with piston from the left
        //col_timesPR =  (xM - xR)/(vR - vM); // time until collision with piston from the right        
        
        for(unsigned i = 0; i < NL; i++){
            double a =  0.5 * g;
            double b = - (vM-vL(i));
            double c = (xL(i)-xM);

            double discrL = b * b - 4 * a * c;

            if(discrL > 0){
                col_timesPLp(i) = (-b + std::pow(discrL, 0.5))/2.0/a;
                col_timesPLm(i) = (-b - std::pow(discrL, 0.5))/2.0/a;
            }
            else{
                col_timesPLp(i) = inf;
                col_timesPLm(i) = inf;
            }
        }

        for(unsigned i = 0; i < NR; i++){
            double a = 0.5 * g;
            double b = - (vM-vR(i));
            double c = (xR(i)-xM);

            double discrR = b * b - 4 * a * c;

            if(discrR > 0){
                col_timesPRp(i) = (-b + std::pow(discrR, 0.5))/2.0/a;
                col_timesPRm(i) = (-b - std::pow(discrR, 0.5))/2.0/a;
            }
            else{
                col_timesPRp(i) = inf;
                col_timesPRm(i) = inf;
            }
        }


#if LOG > 0
        // Update the state vector
        state_t.row(iter)(0) = t;
        state_t.row(iter)(1) = xM;
        state_t.row(iter)(2) = vM;
        for(unsigned i = 0; i < NL; i++){
            state_t.row(iter)(3+i*5) = xL(i);
            state_t.row(iter)(4+i*5) = vL(i);
            state_t.row(iter)(5+i*5) = col_timesL(i);
            state_t.row(iter)(6+i*5) = col_timesPLp(i);
            state_t.row(iter)(7+i*5) = col_timesPLm(i);
        }
        for(unsigned i = NL; i < NL + NR; i++){
            state_t.row(iter)(3+i*5) = xR(i-NL);
            state_t.row(iter)(4+i*5) = vR(i-NL);
            state_t.row(iter)(5+i*5) = col_timesR(i-NL);
            state_t.row(iter)(6+i*5) = col_timesPRp(i-NL);
            state_t.row(iter)(7+i*5) = col_timesPRm(i-NL);
        }
#endif

        // If the particle just colided, it should not be updated, so I simply give it
        // a very large collision time after the calculation
        if(GP_DEBUG) std::cout << "Setting collision time to infinity for particle that just collided\n" << std::flush;
        switch(flag){
            case 0:{col_timesL(pos)   = inf; break; }
            case 1:{col_timesPLm(pos) = inf; break; }
            case 2:{col_timesPLm(pos) = inf; break; }
            case 3:{col_timesR(pos)   = inf; break; }
            case 4:{col_timesPRp(pos) = inf; col_timesPRm(pos) = inf; break; }           
            case 5:{col_timesPRm(pos) = inf; col_timesPRp(pos) = inf; break; }
        }


        // Finding the smallest next collision time
        if(GP_DEBUG) std::cout << "Finding next collision time\n" << std::flush;

        double time_til_col = inf;  // time until the next collision
        double ti;                  // temporary variable

        // Check for collisions with left wall
        for(unsigned i = 0; i < NL; i++){
            ti = col_timesL(i);
            if(ti < time_til_col && ti > 0){
                time_til_col = ti;
                pos = i;
                flag = 0;
            }
        }

        // Check for collisions with piston from the left positive root
        for(unsigned i = 0; i < NL; i++){
            ti = col_timesPLp(i);
            if(ti < time_til_col && ti > 0){
                time_til_col = ti;
                pos = i;
                flag = 1;
            }
        }

        // Check for collisions with piston from the left negative root
        for(unsigned i = 0; i < NL; i++){
            ti = col_timesPLm(i);
            if(ti < time_til_col && ti > 0){
                time_til_col = ti;
                pos = i;
                flag = 2;
            }
        }

        // Check for collisions with right wall
        for(unsigned i = 0; i < NR; i++){
            ti = col_timesR(i);
            if(ti < time_til_col && ti > 0){
                time_til_col = ti;
                pos = i;
                flag = 3;
            }
        }

        // Check for collisions with piston from the right positive root
        for(unsigned i = 0; i < NR; i++){
            ti = col_timesPRp(i);
            if(ti < time_til_col && ti > 0){
                time_til_col = ti;
                pos = i;
                flag = 4;
            }
        }

        // Check for collisions with piston from the right negative root
        for(unsigned i = 0; i < NR; i++){
            ti = col_timesPRm(i);
            if(ti < time_til_col && ti > 0){
                time_til_col = ti;
                pos = i;
                flag = 5;
            }
        }
        

        if(GP_DEBUG) std::cout << t << " " << time_til_col << " " << pos << " " << flag << "\n" << std::flush;


        // Evolve the system in time until the time of the colision
        xL = xL + vL*time_til_col;
        xR = xR + vR*time_til_col;
        xM = xM + vM*time_til_col - 0.5*g*time_til_col*time_til_col;
        vM = vM - g*time_til_col;

        // One of the particles will reverse its velocity
        // This particle should NOT be updated in the next iteration because its col time will be zero
        // Left wall collision
        if(flag == 0){
            vL(pos) *= -1; 
        }

        // collision with piston from the left
        if(flag == 1 || flag ==2){
            double den = 1.0 + gamma;
            double vMtemp = vM;
            double vLtemp = vL(pos);
            vM      = (1.0 - gamma)/den*vMtemp +       2*gamma/den*vLtemp;
            vL(pos) =           2.0/den*vMtemp - (1.0 - gamma)/den*vLtemp;
        }

        // right wall collision
        if(flag == 3){
            vR(pos) *= -1; 
        }

        // collision with piston from the right
        if(flag == 4 || flag == 5){
            double den = 1.0 + gamma;
            double vMtemp = vM;
            double vRtemp = vR(pos);
            vM      = (1.0 - gamma)/den*vMtemp +       2*gamma/den*vRtemp;
            vR(pos) =           2.0/den*vMtemp - (1.0 - gamma)/den*vRtemp;
        }
        t  += time_til_col;
    }

#if LOG > 0
    // Save the log data to a file
    if(GP_DEBUG) std::cout << "Saving to file\n" << std::flush;
    std::ofstream file;
    std::string logname = "log" + std::to_string(identifier) + ".dat";
    file.open(logname);
    file << std::setprecision(9) << state_t;
    file.close();
#endif

    // Save the observables to a file
    //std::ofstream file2;
    //file2.open("observables.dat");
    //file2 << std::setprecision(PRECISION) << observables;
    //file2.close();

    if(GP_DEBUG) std::cout << "Finished\n" << std::flush;
    //return 0;
    return observables;
}
