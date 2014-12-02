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
    timestep.resize(3*this->numberOfBodies());
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
    timestep.resize(3*this->numberOfBodies());
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

std::valarray<double> SolarSystem::calculateForcesAndEnergy(std::valarray<double> X)
{
    forces = 0;
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.setToZero();
    a.setToZero();
    double G = 4*M_PI*M_PI;

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = bodies[j];
            vec3 deltaRVector =  body1.position - body2.position;
            double dr = deltaRVector.length();
            double dx = X[6*i + 0] - X[6*j + 0];
            double dy = X[6*i + 1] - X[6*j + 1];
            double dz = X[6*i + 2] - X[6*j + 2];

            //f is the force acting between body 1 and body 2
            double f = -(G * body1.mass * body2.mass)/(dr*dr*dr);

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

            potentialEnergy -= (body1.mass*body2.mass)/dr;
        }

        kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();

        //angular momentum L = r x v
        double dangularMomentum = body1.position.length() * body1.velocity.length();
        //  std::cout << potentialEnergy << " and KE " << kineticEnergy << "and L" << dangularMomentum<< std::endl;
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

std::valarray<double> SolarSystem::calculateVerlet(std::valarray<double> X, std::valarray<double> A, int bin)
{
    forces = 0;
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.setToZero();
    double G = 4*M_PI*M_PI;

    for(int i=0; i<numberOfBodies(); i++) {
//        std::cout << bodies[i].dtbin << "bin " << bin << std::endl;
        if(bin <= bodies[i].dtbin){
            CelestialBody &body1 = bodies[i];
            for(int j=i+1; j<numberOfBodies(); j++) {
                CelestialBody &body2 = bodies[j];
                vec3 deltaRVector =  body1.position - body2.position;
                double dr = deltaRVector.length();
                double dx = X[3*i + 0] - X[3*j + 0];
                double dy = X[3*i + 1] - X[3*j + 1];
                double dz = X[3*i + 2] - X[3*j + 2];

                //f is the force acting between body 1 and body 2
                double f = -(G * body1.mass * body2.mass)/sqrt(dr*dr*dr);

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

                potentialEnergy -= (body1.mass*body2.mass)/dr;

            }

            kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();

            //angular momentum L = r x v
            //double dangularMomentum = body1.position.length() * body1.velocity.length();
            //std::cout << potentialEnergy << " and KE " << kineticEnergy << "and L" << dangularMomentum<< std::endl;
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

int SolarSystem::numberOfBodies()
{
    return bodies.size();
}

double SolarSystem::totalEnergy()
{
    return this->kineticEnergy + this->potentialEnergy;
}

double SolarSystem::min_time()
{
    // std::valarray<double> timestep(3*numberOfBodies());
    dt_min = .004;
    double acc_constant = 1;
    int bin;
    for(int i =0; i < numberOfBodies();i++){
        double Alength = sqrt(this->A[3*i+0]*this->A[3*i+0] + this->A[3*i+1]*this->A[3*i+1] + this->A[3*i+2]*this->A[3*i+2]);
        double tstep_metric = acc_constant*(1/Alength);

        if (tstep_metric < .01){ bin = 2;}
        else if (tstep_metric <= .5){ bin = 1;}
        else {bin = 0;}
        if(tstep_metric < dt_min && tstep_metric > 0){dt_min = tstep_metric;}
        std::cout << "tstep_metric for body i = " << i << " is " << tstep_metric << " bin = " << bin << std::endl;
        bodies[i].dtbin = bin;
    }
    for(int i = 0; i < numberOfBodies(); i++) {
        CelestialBody &thisBody = bodies[i];
    //std::cout << "bodies.dtbin i= " << i <<  "is " << thisBody.dtbin << std::endl;
    }
   // std::cout << "dt_min is " << dt_min << std::endl;
    return dt_min;
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
std::cout << stepvalues[3*i] << " " << stepvalues[3*i+1]<< " " << bodies[i].dtbin <<std::endl;
    }
    return stepvalues;
}



