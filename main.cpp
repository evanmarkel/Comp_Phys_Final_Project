#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <solarsystem.h>
#include <rk4.h>
#include <verlet.h>
#include <gaussiandeviate.h>
//#include <lib.h>

#include <valarray>
using namespace std;

int main()
{
    ofstream myfile;
    myfile.open("N1000TotalEnergy.txt");

    double G_solarsystem = 4*M_PI*M_PI; //value of Gravitational constant.
    double earthvel = 2*M_PI;

    int N = 1000; //number of particles in simulation
    long seed = -1739562974; //seed number for gaussian random number generator
    double LYtoAU = 63239.73;
    double R = 20; //mean spherical radius to deviate from using random number generator. Units are lightyears
    double mu = 10; //mean particle masses. Relative units are solar masses.

    double rho0 = 3*N*mu/(4*M_PI*R*R*R); //density of the cluster: total mass/volume
    double G = 3*M_PI/(32*rho0); //value of Gravitational constant in relative units. Derived from tau_crunch = 1.
    double tau_crunch = sqrt( (3*M_PI) / (32 * G * rho0) ); //singularity collapse time if continuous fluid. Units of time.

    //eps is the smoothing factor used to reduce the numerical instability of small dr
    //resulting from 2 particles very close together.
    double eps = 0.15;

    cout << "rho" << rho0 << "G" << G << "taucrunch" << tau_crunch << endl;

    //runtime
    double global_min = .005; //step length
    double final_time = 5*tau_crunch; //years

    SolarSystem mySolarSystem;
    for(int i = 0; i < N; i++){
        //random number generator for particle positions
        double sphereu = ran2(&seed);
        double spherev = ran2(&seed);
        double spherew = ran2(&seed);

        double u = sphereu;
        double v = spherev;
        double w = spherew;
        double phi = 2.*M_PI*w;
        double theta = acos(1-2*v);
        double r = R*pow(u,0.333);

        double x = r*sin(theta)*cos(phi);     // transform spherical coordinates into cartesian (x,y,z)
        double y = r*sin(theta)*sin(phi);     // the coordinates represent the positions of the particles. The units are lightyears
        double z = r*cos(theta);              // and the range of values for (x,y,z) are [0,20] lightyears.

        //mass is given with mean 10 Solar Masses with the gaussian deviate returning a normally distributed random number
        //with mean 0 and std deviation of 1
        //small velocities given by random gaussian with standard deviation of .1
        double mean_mass = 10.0; //Solar Masses
        double mass = 10.0 + gaussian_deviate(&seed);

        cout << x << " " << y << " " << z << " " << mass << endl;

        mySolarSystem.addCelestialBody(CelestialBody(vec3(x,y,z), vec3(0,0,0), mass));//vec3(gaussian_deviate(&seed)/1000,gaussian_deviate(&seed)/1000,gaussian_deviate(&seed)/1000), mass));
    }

    //begin runtime calculation
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    //change makeX to makeXV for Verlet. makeX is for the original RK4 valarray with pos, vel combined.
    mySolarSystem.makeXV();
    double step = global_min;

    //timesteps are variable so the system runs a variable number of times, dependent on the steplength each iteration returned by 'step' below
    // and will finish by the 'final_time'
    for(double i = 0; i < final_time; i += step){

        verlet::INTEGRATE(mySolarSystem.X,mySolarSystem.V,mySolarSystem.A,mySolarSystem.E,step*4,mySolarSystem, G, eps);
        step = mySolarSystem.min_time(global_min);
myfile << " total E " << mySolarSystem.CalculateTotalEnergy(mySolarSystem.E) << " i " << i << " step " << step << endl;
        //write to file
       //** myfile << N << endl;
       //** myfile << "stellar cluster in lightyears" << endl;
        for (int i=0;i<mySolarSystem.bodies.size();i++) {
            // myfile << mySolarSystem.X[3*i+0] << " " << mySolarSystem.X[3*i+1] << " " << mySolarSystem.X[3*i+2] << " ";
            //" length" << sqrt(mySolarSystem.A[3*i+0]*mySolarSystem.A[3*i+0] + mySolarSystem.A[3*i+1]*mySolarSystem.A[3*i+1] + mySolarSystem.A[3*i+2]*mySolarSystem.A[3*i+2]) << " ";//Verlet
            //CHANGE X[6*i+0]RK4 to X[3*i+0]Verlet as necessary.
            // myfile << mySolarSystem.X[6*i+0] << " " << mySolarSystem.X[6*i+1] << " " << mySolarSystem.X[6*i+2] << " ";//RK4
            CelestialBody &thisBody = mySolarSystem.bodies[i];
            //myfile << thisBody.mass/60 << " " << mySolarSystem.X[3*i+0] << " " << mySolarSystem.X[3*i+1] << " " << mySolarSystem.X[3*i+2] << endl;
       //**    myfile << "PE " << mySolarSystem.E[3*i] << " KE " << mySolarSystem.E[3*i+1] << " bound " << mySolarSystem.E[3*i+2] << endl;
        }

        //perform RK4 for the timescale of observation
        //change to verlet for verlet
        //RK4::integrate(mySolarSystem.X, global_min, mySolarSystem, G, eps);

    }
    myfile.close();

    //end run time
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    for(int i = 0; i < mySolarSystem.bodies.size(); i++) {
        CelestialBody &thisBody = mySolarSystem.bodies[i];
        cout << "The position of this object is " << thisBody.position << " with velocity " << thisBody.velocity << endl;
    }

    cout << "runtime is " << cpu_time_used << " seconds." << endl;
    return 0;
}

