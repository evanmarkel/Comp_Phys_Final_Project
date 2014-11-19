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
    std::valarray<double> forces;
    double kineticEnergy;
    double potentialEnergy;
    vec3 angularMomentum;
    vec3 a;

    SolarSystem();
    void addCelestialBody(CelestialBody newBody);
    std::valarray<double> calculateForcesAndEnergy(std::valarray<double> X);
    std::valarray<double> calculateVerlet(std::valarray<double> X);
    void makeX();
    void makeXV();
    int numberOfBodies();
    double totalEnergy();
};

#endif // SOLARSYSTEM_H
