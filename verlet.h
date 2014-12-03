#ifndef VERLET_H
#define VERLET_H
#include <solarsystem.h>

#include <vector>
#include <valarray>
class verlet
{
public:
    verlet();
    static void INTEGRATE(std::valarray<double> &X, std::valarray<double> &V, std::valarray<double> &A, double dt, SolarSystem mysystem, double G, double eps);
};

#endif // VERLET_H
