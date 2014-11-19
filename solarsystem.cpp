#include "solarsystem.h"

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

std::valarray<double> SolarSystem::calculateVerlet(std::valarray<double> X)
{
    forces = 0;
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.setToZero();
    double G = 4*M_PI*M_PI;

    for(int i=0; i<numberOfBodies(); i++) {
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
        //double dangularMomentum = body1.position.length() * body1.velocity.length();
        //std::cout << potentialEnergy << " and KE " << kineticEnergy << "and L" << dangularMomentum<< std::endl;
    }
    for (int i =0; i <numberOfBodies(); i++){
        double m = bodies[i].mass;
        A[3*i + 0] = forces[3*i + 0] / m; //vx(t+dt/2)=vx(t)+.5ax(t)dt
        A[3*i + 1] = forces[3*i + 1] / m; //vy(t+dt/2)=vy(t)+.5ay(t)dt
        A[3*i + 2] = forces[3*i + 2] / m; //vz(t+dt/2)=vz(t)+.5az(t)dt
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