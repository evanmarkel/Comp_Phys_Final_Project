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
    myfile.open("es.txt");
    //define step variable h and run length. years are the units
    double h = .004; //step length
    double final_time = 400; //years
    double G = 4*M_PI*M_PI; //value of Gravitational constant.
    double earthvel = 2*M_PI;

    int N = 1; //number of particles in simulation
    long seed = -1739562974; //seed number for gaussian random number generator
    double LYtoAU = 63239.73;
    double R = 20; //mean spherical radius to deviate from using random number generator. Units are lightyears
    
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
    double y = r*sin(theta)*sin(phi);
    double z = r*cos(theta);

    //mass is given with mean 10 Solar Masses with the gaussian deviate returning a normally distributed random number
    //with mean 0 and std deviation of 1
    //small velocities given by random gaussian with standard deviation of .1
    double mean_mass = 10.0; //Solar Masses
    double mass = 10.0 + gaussian_deviate(&seed);

    cout << x << " " << y << " " << z << " " << mass << endl;

      //  mySolarSystem.addCelestialBody(CelestialBody(vec3(x,y,z), vec3(gaussian_deviate(&seed)/10,gaussian_deviate(&seed)/10,gaussian_deviate(&seed)/10), mass));
    }
    CelestialBody sun(vec3(0,0,0), vec3(0,0,0), 1.0);
    //CelestialBody mercury(vec3(-.39,0,0), vec3(0,-1.67*earthvel,0), 1.65e-7);
    //CelestialBody venus(vec3(0,.72,0), vec3(1.174*earthvel,0,0), 2.44e-6);
    CelestialBody earth(vec3(0,1,0), vec3(earthvel,0, 0), 3e-6);
    //CelestialBody earthmoon(vec3(1.00256956, 0, 0), vec3(0, .034352+earthvel, 0), 3.6e-8);
    //CelestialBody mars(vec3(0,1.52,0), vec3(.802*earthvel,0, 0), 0.0000003209425513);
    //CelestialBody jupiter(vec3(5.2, 0, 0), vec3(0, .434*earthvel, 0), 3e-3);  //Jupiter's orbital velocity is 43% of earth's 1/r^2. or 5.20/12year orbit.
    //CelestialBody saturn(vec3(0,-9.54,0), vec3(-.323*earthvel,0, 0), 2.85e-5);
    //CelestialBody uranus(vec3(19.19,0,0), vec3(0,.228*earthvel, 0), 4.3e-5);
    //CelestialBody neptune(vec3(30.06,0,0), vec3(0,.182*earthvel, 0), 5.516e-5);
    //CelestialBody pluto(vec3(0,39.53,0), vec3(.159*earthvel,0, 0), 7.79e-9);


    mySolarSystem.addCelestialBody(sun);
    //mySolarSystem.addCelestialBody(mercury);
    //mySolarSystem.addCelestialBody(venus);
    mySolarSystem.addCelestialBody(earth);
    //mySolarSystem.addCelestialBody(earthmoon);
    //mySolarSystem.addCelestialBody(mars);
    //mySolarSystem.addCelestialBody(jupiter);
    //mySolarSystem.addCelestialBody(saturn);
    //mySolarSystem.addCelestialBody(uranus);
    //mySolarSystem.addCelestialBody(neptune);
    //mySolarSystem.addCelestialBody(pluto);

    //begin runtime check
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    //change makeX to makeXV for Verlet. makeX is for the original RK4 valarray with pos, vel combined.
    mySolarSystem.makeXV();
    double step = h;
    for(int i = 0; i < final_time/h; i++){

        verlet::INTEGRATE(mySolarSystem.X,mySolarSystem.V,mySolarSystem.A,step,mySolarSystem);
        step = mySolarSystem.min_time();

        //write to file
        //CHANGE X[6*i+0]RK4 to X[3*i+0]Verlet
        for (int i=0;i<mySolarSystem.bodies.size();i++) {
            myfile << mySolarSystem.X[3*i+0] << " " << mySolarSystem.X[3*i+1] << " " << mySolarSystem.X[3*i+2] << " ";//" length" << sqrt(mySolarSystem.A[3*i+0]*mySolarSystem.A[3*i+0] + mySolarSystem.A[3*i+1]*mySolarSystem.A[3*i+1] + mySolarSystem.A[3*i+2]*mySolarSystem.A[3*i+2]) << " ";//Verlet
           // myfile << mySolarSystem.X[6*i+0] << " " << mySolarSystem.X[6*i+1] << " " << mySolarSystem.X[6*i+2] << " ";//RK4
        }
        myfile << endl;

        //perform RK4 for the timescale of observation
        //change to verlet for verlet
        //RK4::integrate(mySolarSystem.X,h,mySolarSystem);



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

