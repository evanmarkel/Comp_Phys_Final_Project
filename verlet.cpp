#include "verlet.h"

verlet::verlet()
{
}

void verlet::INTEGRATE(std::valarray<double> &X, std::valarray<double> &V, std::valarray<double> &A , double dt, SolarSystem mysystem, double G, double eps)
{

    for(int i=0; i<4; i++) {
        if(i==0) {
            V = V + 0.5*dt*A*mysystem.bin_particles(0); //1
            V = V + 0.5*dt*0.5*A*mysystem.bin_particles(1); //1
            V = V + 0.5*dt*0.25*A*mysystem.bin_particles(2); //1

            X = X + V*dt*0.25; //2
            A = mysystem.calculateVerlet(X,A,2,G,eps); //3
            V = V + 0.5*dt*0.25*A*mysystem.bin_particles(2); //4
        } else if(i==1) {
            V = V + 0.5*dt*0.25*A*mysystem.bin_particles(2); //1

            X = X + V*dt*0.25; //2
            A = mysystem.calculateVerlet(X,A,1,G,eps); //3
            V = V + 0.5*dt*0.5*A*mysystem.bin_particles(1); //4
            V = V + 0.5*dt*0.25*A*mysystem.bin_particles(2); //4
        } else if(i==2) {
            V = V + 0.5*dt*0.5*A*mysystem.bin_particles(1); //1
            V = V + 0.5*dt*0.25*A*mysystem.bin_particles(2); //1

            X = X + V*dt*0.25; //2
            A = mysystem.calculateVerlet(X,A,2,G,eps); //3
            V = V + 0.5*dt*0.25*A*mysystem.bin_particles(2); //4
        } else {
            V = V + 0.5*dt*0.25*A*mysystem.bin_particles(2); //1

            X = X + V*dt*0.25; //2
            A = mysystem.calculateVerlet(X,A,0,G,eps); //3
            V = V + 0.5*A*dt*mysystem.bin_particles(0); //4
            V = V + 0.5*0.5*A*dt*mysystem.bin_particles(1); //4
            V = V + 0.5*0.25*A*dt*mysystem.bin_particles(2); //4

        }

    }

}
