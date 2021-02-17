#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include "gravity_piston.hpp"

#define inf 999999.9

// Compilar com 
// g++ -O3 -Wall -fopenmp -I/usr/include/eigen3 -DLOG=0 gravity_piston.cpp -o piston

Eigen::Array<double, -1, -1> gravity_piston(unsigned NL, unsigned NR, unsigned NIter, double gamma, double betat_L, double betat_R, unsigned identifier){
    // Arguments to this function:
    //         NL - number of particles in the left chamber of the piston
    //         NR - number of particles in the right chamber of the piston
    //      NIter - number of iterations we want, or number of collisions
    //      gamma - mass of one particle divided by the mass of the wall
    //         TL - initial temperature of the left chamber
    //         TR - initial temperature of the right chamber
    // identifier - unique identifier assigned by the paralelization

    Eigen::Array<double, -1, -1> xL(NL,1), xR(NR,1), vL(NL,1), vR(NR,1); // Positions and velocities of the particles
    Eigen::Array<double, -1, -1> col_timesL(NL,1), col_timesR(NR,1);     // Collision times with left and right walls
    Eigen::Array<double, -1, -1> col_timesPLp(NL,1), col_timesPLm(NL,1), col_timesPRp(NR,1), col_timesPRm(NR,1); // collision times with piston

    double xM, vM; // Instantaneous position and velocity of the piston
    Eigen::Array<double, -1, -1> observables(NIter, 5); // Tracking some variables over time, such as xM


    // List of times, positions, velocities and collision times for every particle.
    // This generates a VERY large file, so it is most useful for debugging and
    // is therefore an optional compile-time flag
#if LOG > 0
    Eigen::Array<double, -1, -1> state_t;
    state_t = Eigen::Array<double, -1, -1>::Zero(NIter, 5*(NL + NR) + 3); // State vector for all iterations
#endif


    // Initializing positions
    xM = 0.5; // piston starts in the middle
    xL = (Eigen::Array<double, -1, -1>::Random(NL, 1) + 1.0)/2.0*xM*0.99;
    xR = (Eigen::Array<double, -1, -1>::Random(NR, 1) + 1.0)/2.0*(1-xM)*0.99 + xM;

    // Initializing the velocities uniformly
    vM = 0.; // piston starts without velocity
    vL = Eigen::Array<double, -1, -1>::Random(NL, 1);
    vR = Eigen::Array<double, -1, -1>::Random(NR, 1);

    // Initializing the velocities according to the Maxwell-Boltzmann distribution
    std::mt19937 rnd;
    std::normal_distribution<double> normal_dist{0., 1. };
    std::random_device r;
    std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()};
    rnd.seed(seed2);

    for(unsigned i = 0; i < NL; i++){
      	vL(i) = normal_dist(rnd)/std::sqrt(2*betat_L*gamma);
    } 

    for(unsigned i = 0; i < NR; i++){
      	vR(i) = normal_dist(rnd)/std::sqrt(2*betat_R*gamma);
    } 

    double vM_temp;
    vM_temp = vM; // velocity of the piston before the collision

    double t = 0; // time variable

    int flag = -1;     // 0: left wall col, 1: left piston col, 2: right wall col, 3: right pison col
    unsigned pos = -1; // index of the particle that will collide. can be a particle on the left or right

    for(unsigned iter = 0; iter < NIter; iter++){
        if(GP_DEBUG) std::cout << "ITERATION: " << iter << "\n" << std::flush;

        observables(iter,0) = t;
        observables(iter,1) = xM;
        observables(iter,2) = (vL * vL).mean(); // Measure of the average kinetic energy of the left gas
        observables(iter,3) = (vR * vR).mean(); // Measure of the average kinetic energy of the right gas
        observables(iter,4) = gamma * (vM-vM_temp); // Measure of the momentum transfered 


        // calculate all the times until collision with left and right walls
        col_timesL  = -(xL - 0.0)/vL; // time until collision with left wall
        col_timesR  =  (1.0 - xR)/vR; // time until collision with right wall
        
        // Calculate all the collision times with piston from the left
        // This requires solving a quadratic equation which may have no solution. In that case
        // the time until collision is set to infinity
        for(unsigned i = 0; i < NL; i++){
            double a =  0.25;
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

        // Calculate all the collision times with piston from the right
        for(unsigned i = 0; i < NR; i++){
            double a = 0.25;
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

        // If the particle just colided, special care has to be taken, because in the next
        // iteration, it will be in exactly the same site as the object with which it collided,
        // so the time for the next collision with it will be zero. This next collision therefore
        // has to be checked manually and discarded appropriately.
        //
        // If it collided with the wall, it cannot collide with it again in the next iteration. It will either 
        // collide with the piston or another particle collides somewhere. So, if it collided with the 
        // wall, the collision time for this (which is zero) is set manually to infinity.
        //
        // If it collided with the piston from the left, since the piston has gravity, it may collide
        // with it AGAIN in the next iteration. In this circumstance, since there will be two possible times
        // for the collision (parabola intersecting line), the earliest one will necessarily be zero and the
        // second one not, so the earliest collision is set to infinity
        //
        // If it collided with the piston from the right, it simply cannot collide again in the next iteration.
        if(GP_DEBUG) std::cout << "Setting collision time to infinity for particle that just collided\n" << std::flush;
        switch(flag){
            case 0:{col_timesL(pos)   = inf; break; }
            case 1:{col_timesPLm(pos) = inf; break; }
            case 2:{col_timesPLm(pos) = inf; break; }
            case 3:{col_timesR(pos)   = inf; break; }
            case 4:{col_timesPRp(pos) = inf; col_timesPRm(pos) = inf; break; }           
            case 5:{col_timesPRm(pos) = inf; col_timesPRp(pos) = inf; break; }
        }


        // Finding the smallest next collision time (time_til_col), the particle that will collide (pos) and
        // the type of collision that is about to happen (flag). 
        //    flag = 0: Left particle will collide with left wall
        //    flag = 1: Left particle will collide with piston (positive root)
        //    flag = 2: Left particle will collide with piston (negative root)
        //    flag = 3: Right particle will collide with right wall
        //    flag = 4: Right particle will collide with piston (positive root)
        //    flag = 5: Right particle will collide with piston (negative root)
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
        xM = xM + vM*time_til_col - 0.25*time_til_col*time_til_col;
        vM = vM - 0.5*time_til_col;

        vM_temp = vM; // Velocity of the piston before the collision

        // Perform the collision according to the type of collision that will happen
        // The velocities have to be updated

        // Left wall collision
        if(flag == 0){
            vL(pos) *= -1; 
        }

        // Elastic collision with piston from the left
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

        // Elastic collision with piston from the right
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

    if(GP_DEBUG) std::cout << "Finished\n" << std::flush;
    return observables;
}
