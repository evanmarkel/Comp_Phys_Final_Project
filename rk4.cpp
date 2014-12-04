#include "rk4.h"

RK4::RK4()
{
}

void RK4::integrate(std::valarray<double> &X, std::valarray<double> &V, double dt, SolarSystem mysystem, double G, double eps)
{
    std::valarray<double> k1(1,6*mysystem.numberOfBodies());
    std::valarray<double> k2(1,6*mysystem.numberOfBodies());
    std::valarray<double> k3(1,6*mysystem.numberOfBodies());
    std::valarray<double> k4(1,6*mysystem.numberOfBodies());

    // RK4 integration using vector X from solarysystem class.
    k1 = mysystem.calculateRK4(X, V, G, eps) * dt;
    k2 = mysystem.calculateRK4(X + 0.5 * k1, V, G, eps) * dt;
    k3 = mysystem.calculateRK4(X + 0.5 * k2, V, G, eps) * dt;
    k4 = mysystem.calculateRK4(X + k3, V, G, eps) * dt;
    X += (1.0/6) * (k1 + 2 * (k2 + k3) + k4);

}
