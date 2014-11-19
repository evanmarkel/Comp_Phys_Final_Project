#include "rk4.h"

RK4::RK4()
{
}

void RK4::integrate(std::valarray<double> &X, double dt, SolarSystem mysystem)
{
    std::valarray<double> k1(1,6*mysystem.numberOfBodies());
    std::valarray<double> k2(1,6*mysystem.numberOfBodies());
    std::valarray<double> k3(1,6*mysystem.numberOfBodies());
    std::valarray<double> k4(1,6*mysystem.numberOfBodies());

    // RK4 integration using vector X from solarysystem class.
    k1 = mysystem.calculateForcesAndEnergy(X) * dt;
    k2 = mysystem.calculateForcesAndEnergy(X + 0.5 * k1) * dt;
    k3 = mysystem.calculateForcesAndEnergy(X + 0.5 * k2) * dt;
    k4 = mysystem.calculateForcesAndEnergy(X + k3) * dt;
    X += (1.0/6) * (k1 + 2 * (k2 + k3) + k4);

}
