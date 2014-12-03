#ifndef RK4_H
#define RK4_H
#include <solarsystem.h>

#include <vector>
#include <valarray>
using std::vector;

class RK4
{
public:
    RK4();
    static void integrate(std::valarray<double> &X, double dt, SolarSystem mysystem, double G, double eps);

};

#endif // RK4_H
