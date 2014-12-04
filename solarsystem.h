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
    std::valarray<double> E;
    std::valarray<double> A;
    std::valarray<double> V;
    std::valarray<double> k;
    std::valarray<double> v;
    std::valarray<double> forces;
    std::valarray<double> stepvalues;
    double kineticEnergy;
    double potentialEnergy;
    double timeclass;
    double dt_min;
    double array_time;
    int isBound;
    vec3 angularMomentum;
    vec3 a;

    SolarSystem();
    void addCelestialBody(CelestialBody newBody);
    std::valarray<double> calculateRK4(std::valarray<double> X, std::valarray<double> V, double G, double eps);
    std::valarray<double> calculateVerlet(std::valarray<double> X,std::valarray<double> A, std::valarray<double> V, int bin, double G, double eps);
    std::valarray<double> calculateEnergy(std::valarray<double> X, std::valarray<double> V, double G);
    void makeX();
    void makeXV();
    int numberOfBodies();
    double CalculateTotalEnergy(std::valarray<double> E);
    double min_time(double global_min);
    std::valarray<double> bin_particles(int boolean_bin);
};

#endif // SOLARSYSTEM_H
