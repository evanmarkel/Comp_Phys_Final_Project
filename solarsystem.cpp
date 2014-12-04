#include "solarsystem.h"
#include <math.h>

SolarSystem::SolarSystem()
{
}

void SolarSystem::addCelestialBody(CelestialBody newBody) {
    bodies.push_back(newBody);
}

void SolarSystem::makeX(){

    X.resize(6*this->numberOfBodies());
    forces.resize(3*this->numberOfBodies());
    k.resize(6*this->numberOfBodies());
    v.resize(6*this->numberOfBodies());
    E.resize(3*this->numberOfBodies());
    for(int i = 0; i<this->numberOfBodies(); i++){
        this->X[6*i + 0] = this->bodies[i].position[0];
        this->X[6*i + 1] = this->bodies[i].position[1];
        this->X[6*i + 2] = this->bodies[i].position[2];
        this->X[6*i + 3] = this->bodies[i].velocity[0];
        this->X[6*i + 4] = this->bodies[i].velocity[1];
        this->X[6*i + 5] = this->bodies[i].velocity[2];
    }
}

void SolarSystem::makeXV(){

    X.resize(3*this->numberOfBodies());
    A.resize(3*this->numberOfBodies());
    V.resize(3*this->numberOfBodies());
    E.resize(3*this->numberOfBodies());
    forces.resize(3*this->numberOfBodies());
    for(int i = 0; i<this->numberOfBodies(); i++){
        this->X[3*i + 0] = this->bodies[i].position[0];
        this->X[3*i + 1] = this->bodies[i].position[1];
        this->X[3*i + 2] = this->bodies[i].position[2];
        this->V[3*i + 0] = this->bodies[i].velocity[0];
        this->V[3*i + 1] = this->bodies[i].velocity[1];
        this->V[3*i + 2] = this->bodies[i].velocity[2];
    }
}

std::valarray<double> SolarSystem::calculateRK4(std::valarray<double> X, std::valarray<double> V, double G, double eps)
{
    forces = 0;
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.setToZero();
    a.setToZero();

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = bodies[j];

            double dx = X[6*i + 0] - X[6*j + 0];
            double dy = X[6*i + 1] - X[6*j + 1];
            double dz = X[6*i + 2] - X[6*j + 2];

            double dr = sqrt((dx*dx + dy*dy + dz*dz));

            //f is the force acting between body 1 and body 2
            double f = -(G * body1.mass * body2.mass)/((dr*dr+ eps*eps)*dr);

            //a is the force multiplied by the relative position of the two bodies. the movement is added to previous position
            double axtemp = dx*(f);
            double aytemp = dy*(f);
            double aztemp = dz*(f);

            //forces vector values for body 1
            forces[3*i + 0] += axtemp;
            forces[3*i + 1] += aytemp;
            forces[3*i + 2] += aztemp;

            //forces vector values for body 2
            forces[3*j + 0] -= axtemp;
            forces[3*j + 1] -= aytemp;
            forces[3*j + 2] -= aztemp;

            //potentialEnergy -= (body1.mass*body2.mass)/dr;
        }

        //kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();

    }
    //calculate Vx Vy Vz.assign to X. also assign forces.
    for (int i =0; i <numberOfBodies(); i++){
        double m = bodies[i].mass;
        k[6*i + 0] = X[6*i +3]; //vx
        k[6*i + 1] = X[6*i +4]; //vy
        k[6*i + 2] = X[6*i +5]; //vy
        k[6*i + 3] = forces[3*i + 0] / m; //Fx/m
        k[6*i + 4] = forces[3*i + 1] / m; //Fy/m
        k[6*i + 5] = forces[3*i + 2] / m; //Fz/m
    }
    return k;
}

std::valarray<double> SolarSystem::calculateVerlet(std::valarray<double> X, std::valarray<double> A, std::valarray<double> V, int bin, double G, double eps)
{
    forces = 0;
    isBound = 2;
    angularMomentum.setToZero();

    for(int i=0; i<numberOfBodies(); i++) {
        kineticEnergy = 0;
        potentialEnergy = 0;
        if(bin <= bodies[i].dtbin){
            CelestialBody &body1 = bodies[i];
            for(int j=i+1; j<numberOfBodies(); j++) {
                CelestialBody &body2 = bodies[j];

                // parallelize
                //#pragma omp parallel for shared(dx,dy,dz) private(dr)
                double dx = X[3*i + 0] - X[3*j + 0];
                double dy = X[3*i + 1] - X[3*j + 1];
                double dz = X[3*i + 2] - X[3*j + 2];

                double dr = sqrt((dx*dx + dy*dy + dz*dz));

                //f is the force acting between body 1 and body 2
                double f = -(G * body1.mass * body2.mass)/((dr*dr+ eps*eps)*dr);

                double axtemp = dx*(f);
                double aytemp = dy*(f);
                double aztemp = dz*(f);

                //forces vector values for body 1
                forces[3*i + 0] += axtemp;
                forces[3*i + 1] += aytemp;
                forces[3*i + 2] += aztemp;

                //forces vector values for body 2
                forces[3*j + 0] -= axtemp;
                forces[3*j + 1] -= aytemp;
                forces[3*j + 2] -= aztemp;
            }
        }
    }

    for (int i =0; i <numberOfBodies(); i++){
        if(bin <= bodies[i].dtbin){
            double m = bodies[i].mass;
            A[3*i + 0] = forces[3*i + 0] / m;
            A[3*i + 1] = forces[3*i + 1] / m;
            A[3*i + 2] = forces[3*i + 2] / m;
        }
    }
    return A;
}

std::valarray<double> SolarSystem::calculateEnergy(std::valarray<double> X, std::valarray<double> V, double G)
{
    isBound = 2;
    angularMomentum.setToZero();

    for(int i=0; i<numberOfBodies(); i++) {
        kineticEnergy = 0;
        potentialEnergy = 0;

        CelestialBody &body1 = bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = bodies[j];

            double dx = X[3*i + 0] - X[3*j + 0];
            double dy = X[3*i + 1] - X[3*j + 1];
            double dz = X[3*i + 2] - X[3*j + 2];

            double dr = sqrt((dx*dx + dy*dy + dz*dz));

            potentialEnergy -= G * (body1.mass*body2.mass) / (dr);
        }
        //calcuate v^2 for Kinetic Energy calculation.
        double dvx = V[3*i + 0];
        double dvy = V[3*i + 1];
        double dvz = V[3*i + 2];

        double dv = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

        //angular momentum L = r x v
        double dx1 = X[3*i + 0];
        double dy1 = X[3*i + 1];
        double dz1 = X[3*i + 2];

        double dr1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        double dangularMomentum = dr1 * dv;

        //KE = (1/2) * Mass * Velocity^2
        kineticEnergy = 0.5*body1.mass*dv*dv*.72;

        E[3*i + 0] = potentialEnergy;
        E[3*i + 1] = kineticEnergy;

        //test if 2KE=U holds and satisfies virial theorem. Particle can escape and is ejected if K - U > 0 and otherwise bound.
        if (potentialEnergy + kineticEnergy < 0) { isBound = 1;}
        else { isBound = 0;}
        E[3*i+2] = isBound;
    }
    return E;
}

int SolarSystem::numberOfBodies()
{
    return bodies.size();
}

double SolarSystem::CalculateTotalEnergy(std::valarray<double> E)
{
    //calculates the total energy for the bound particles
    double totalkinetic = 0;
    double totalpotential = 0;
    for(int i=0; i<numberOfBodies(); i++) {
        totalpotential += E[3*i]*E[3*i + 2]/2;
        totalkinetic += E[3*i + 1]*E[3*i + 2];
    }
    std::cout << "virial check " << 2*totalkinetic/totalpotential << std::endl;
    return totalkinetic + totalpotential;
}

double SolarSystem::min_time(double global_min)
{
    dt_min = global_min*4;
    double acc_constant = 1;
    int bin;
    double temp_step = dt_min;
    for(int i =0; i < numberOfBodies();i++){
        double Alength = sqrt(this->A[3*i+0]*this->A[3*i+0] + this->A[3*i+1]*this->A[3*i+1] + this->A[3*i+2]*this->A[3*i+2]);
        double tstep_metric = acc_constant*(1/Alength);

        if (tstep_metric < .0012){ bin = 2; dt_min = global_min;}
        else if (tstep_metric <= .1){ bin = 1; dt_min = global_min*2;}
        else {bin = 0; dt_min = global_min*4;}
        // if(tstep_metric < dt_min && tstep_metric > global_min){dt_min = dt_min;}
        std::cout << "tstep_metric for body i = " << i << " is " << tstep_metric << " bin = " << bin << std::endl;
        bodies[i].dtbin = bin;
        if(dt_min < temp_step){ temp_step = dt_min;}
    }
    for(int i = 0; i < numberOfBodies(); i++) {
        CelestialBody &thisBody = bodies[i];
        //std::cout << "bodies.dtbin i= " << i <<  "is " << thisBody.dtbin << std::endl;
    }
    std::cout << "dt_min is " << temp_step << std::endl;
    return temp_step;
}

std::valarray<double> SolarSystem::bin_particles(int boolean_bin)
{
    stepvalues.resize(3*this->numberOfBodies());
    for(int i=0;i < numberOfBodies(); i++){
        if (bodies[i].dtbin == boolean_bin){
            stepvalues[3*i + 0]=1;
            stepvalues[3*i + 1]=1;
            stepvalues[3*i + 2]=1;
        }
        else{
            stepvalues[3*i + 0]=0;
            stepvalues[3*i + 1]=0;
            stepvalues[3*i + 2]=0;
        }
    }
    return stepvalues;
}



