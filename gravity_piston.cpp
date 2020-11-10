#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <iomanip>
#include <random>

#define inf 999999.9
#define DEBUG 0
#define PRECISION 9
#define KB 1
//#define LOG 0

// Verificar edge cases desta implementação
//
// Compilar com 
// g++ -Wall -fopenmp -I/usr/include/eigen3 -DLOG=0 gravity_piston.cpp -o piston

int main(int argc, char** argv){

    // Get the command line parameters
    unsigned NL, NR, NIter;
    double gamma, TL, TR;
    if(argc != 7){
        std::cout << "Wrong number of parameters. Exiting.\n";
        exit(1);
    }

    NL    = atoi(argv[1]);
    NR    = atoi(argv[2]);
    NIter = atoi(argv[3]);
    gamma = atof(argv[4]);
    TL    = atof(argv[5]);
    TR    = atof(argv[6]);

    std::cout << "Program parameters:\n";
    std::cout << "NL: " << NL << "\n";
    std::cout << "NR: " << NR << "\n";
    std::cout << "NIter: " << NIter << "\n";
    std::cout << "gamma: " << gamma << "\n";

    Eigen::Array<double, -1, -1> xL(NL,1), xR(NR,1), vL(NL,1), vR(NR,1); //positions and velocities
    Eigen::Array<double, -1, -1> col_timesL(NL,1), col_timesR(NR,1), col_timesPL(NL,1), col_timesPR(NR,1); // collision times
    double xM, vM; // Instantaneous position and velocity of the piston
    Eigen::Array<double, -1, -1> observables(NIter, 5); // tracking some variables across time suhc as xM
    double vM_temp;


    // list of times, positions, velocities and collision times for every particle
#if LOG > 0
    Eigen::Array<double, -1, -1> state_t;
    state_t = Eigen::Array<double, -1, -1>::Zero(NIter, 4*(NL + NR) + 3); // State vector for all iterations
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
        if(DEBUG) std::cout << "ITERATION: " << iter << "\n" << std::flush;

        observables(iter,0) = t;
        observables(iter,1) = xM;
        observables(iter,2) = (vL * vL).mean();
        observables(iter,3) = (vR * vR).mean();
        observables(iter,4) = gamma * (vM-vM_temp);

        vM_temp = vM;

        // calculate all the colision times
        col_timesL  = -(xL - 0.0)/vL;       // time until collision with left wall
        col_timesR  =  (1.0 - xR)/vR;       // time until collision with right wall
        col_timesPL =  (xM - xL)/(vL - vM); // time until collision with piston from the left
        col_timesPR =  (xM - xR)/(vR - vM); // time until collision with piston from the right


#if LOG > 0
        // Update the state vector
        state_t.row(iter)(0) = t;
        state_t.row(iter)(1) = xM;
        state_t.row(iter)(2) = vM;
        for(unsigned i = 0; i < NL; i++){
            state_t.row(iter)(3+i*4) = xL(i);
            state_t.row(iter)(4+i*4) = vL(i);
            state_t.row(iter)(5+i*4) = col_timesL(i);
            state_t.row(iter)(6+i*4) = col_timesPL(i);
        }
        for(unsigned i = NL; i < NL + NR; i++){
            state_t.row(iter)(3+i*4) = xR(i-NL);
            state_t.row(iter)(4+i*4) = vR(i-NL);
            state_t.row(iter)(5+i*4) = col_timesR(i-NL);
            state_t.row(iter)(6+i*4) = col_timesPR(i-NL);
        }
#endif

        // If the particle just colided, it should not be updated, so I simply give it
        // a very large collision time after the calculation
        if(DEBUG) std::cout << "Setting collision time to infinity for particle that just collided\n" << std::flush;
        switch(flag){
            case 0:{col_timesL(pos)  = inf; break; }
            case 1:{col_timesPL(pos) = inf; break; }
            case 2:{col_timesR(pos)  = inf; break; }
            case 3:{col_timesPR(pos) = inf; break; }
        }


        // Finding the smallest next collision time
        if(DEBUG) std::cout << "Finding next collision time\n" << std::flush;

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

        // Check for collisions with piston from the left
        for(unsigned i = 0; i < NL; i++){
            ti = col_timesPL(i);
            if(ti < time_til_col && ti > 0){
                time_til_col = ti;
                pos = i;
                flag = 1;
            }
        }

        // Check for collisions with right wall
        for(unsigned i = 0; i < NR; i++){
            ti = col_timesR(i);
            if(ti < time_til_col && ti > 0){
                time_til_col = ti;
                pos = i;
                flag = 2;
            }
        }

        // Check for collisions with piston from the right
        for(unsigned i = 0; i < NR; i++){
            ti = col_timesPR(i);
            if(ti < time_til_col && ti > 0){
                time_til_col = ti;
                pos = i;
                flag = 3;
            }
        }

        if(DEBUG) std::cout << t << " " << time_til_col << " " << pos << " " << flag << "\n" << std::flush;


        // Evolve the system in time until the time of the colision
        xL = xL + vL*time_til_col;
        xR = xR + vR*time_til_col;
        xM = xM + vM*time_til_col;

        // One of the particles will reverse its velocity
        // This particle should NOT be updated in the next iteration because its col time will be zero
        // Left wall collision
        if(flag == 0){
            vL(pos) *= -1; 
        }

        // collision with piston from the left
        if(flag == 1){
            double den = 1.0 + gamma;
            double vMtemp = vM;
            double vLtemp = vL(pos);
            vM      = (1.0 - gamma)/den*vMtemp +       2*gamma/den*vLtemp;
            vL(pos) =           2.0/den*vMtemp - (1.0 - gamma)/den*vLtemp;
        }

        // right wall collision
        if(flag == 2){
            vR(pos) *= -1; 
        }

        // collision with piston from the right
        if(flag == 3){
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
    if(DEBUG) std::cout << "Saving to file\n" << std::flush;
    std::ofstream file;
    file.open("log.dat");
    file << std::setprecision(PRECISION) << state_t;
    file.close();
#endif

    // Save the observables to a file
    std::ofstream file2;
    file2.open("observables.dat");
    file2 << std::setprecision(PRECISION) << observables;
    file2.close();

    if(DEBUG) std::cout << "Finished\n" << std::flush;
    return 0;
}
