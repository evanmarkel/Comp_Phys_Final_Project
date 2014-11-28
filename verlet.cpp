#include "verlet.h"

verlet::verlet()
{
}

void verlet::INTEGRATE(std::valarray<double> &X, std::valarray<double> &V, std::valarray<double> &A , double dt, SolarSystem mysystem)
{

    for(int i=0; i<4; i++) {
        if(i==0) {
            V = V + 0.5*mysystem.dt_min*4*A*mysystem.bin_particles(0);
            V = V + 0.5*mysystem.dt_min*2*A*mysystem.bin_particles(1);
            V = V + 0.5*mysystem.dt_min*A*mysystem.bin_particles(2);

            X = X + V*mysystem.dt_min;
            A = mysystem.calculateVerlet(X,A,2);
            V = V + 0.5*mysystem.dt_min*A*mysystem.bin_particles(2);
        } else if(i==1) {
            V = V + 0.5*mysystem.dt_min*2*A*mysystem.bin_particles(1);
            V = V + 0.5*mysystem.dt_min*A*mysystem.bin_particles(2);

            X = X + V*mysystem.dt_min;
            A = mysystem.calculateVerlet(X,A,1);
            V = V + 0.5*mysystem.dt_min*A*mysystem.bin_particles(1);
        } else if(i==2) {
            V = V + 0.5*mysystem.dt_min*A*mysystem.bin_particles(2);

            X = X + V*mysystem.dt_min;
            A = mysystem.calculateVerlet(X,A,2);
            V = V + 0.5*mysystem.dt_min*A*mysystem.bin_particles(2);
        } else {
            V = V + 0.5*mysystem.dt_min*4*A*mysystem.bin_particles(0);
            V = V + 0.5*mysystem.dt_min*2*A*mysystem.bin_particles(1);
            V = V + 0.5*mysystem.dt_min*A*mysystem.bin_particles(2);

            X = X + V*mysystem.dt_min;
            A = mysystem.calculateVerlet(X,A,0);
            V = V + 0.5*A*mysystem.dt_min*mysystem.bin_particles(0);

        }

    }

}
