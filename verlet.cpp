#include "verlet.h"

verlet::verlet()
{
}

void verlet::INTEGRATE(std::valarray<double> &X, std::valarray<double> &V, std::valarray<double> &A , double dt, SolarSystem mysystem)
{
    V = V + 0.5*dt *A;
    X = X + V*dt;
    A = mysystem.calculateVerlet(X);
    V = V + 0.5*A*dt;
}
