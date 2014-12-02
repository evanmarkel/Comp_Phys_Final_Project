#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include <celestialbody.h>
#include <vector>
#include <valarray>

using std::vector;

class SolarSystem
{
public:
    vector<CelestialBody> bodies;
    std::valarray<double> X;
    std::valarray<double> A;
    std::valarray<double> V;
    std::valarray<double> k;
    std::valarray<double> v;
    std::valarray<double> timestep;
    std::valarray<double> forces;
    std::valarray<double> stepvalues;
    double kineticEnergy;
    double potentialEnergy;
    double timeclass;
    double dt_min;
    double array_time;
    vec3 angularMomentum;
    vec3 a;

    SolarSystem();
    void addCelestialBody(CelestialBody newBody);
    std::valarray<double> calculateForcesAndEnergy(std::valarray<double> X);
    std::valarray<double> calculateVerlet(std::valarray<double> X,std::valarray<double> A, int bin);
    void makeX();
    void makeXV();
    int numberOfBodies();
    double totalEnergy();
    double min_time();
    std::valarray<double> bin_particles(int boolean_bin);
};

#endif // SOLARSYSTEM_H
