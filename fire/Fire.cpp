//
// Created by Matthew Nicoletti on 2019-05-07.
//

#include "Fire.h"







void Fire::update() {

    // solve advection



}

void Fire::propagateFront(double w1, double w2, double w3) {
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            for (int k = 1; k <= N; k++) {
                double nx, ny, nz;
                // components of the gradient
                nx = (grid[4*((i + 1)*N*N + j*N       + k)]     - grid[4*((i - 1)*N*N + j*N       + k)]) / (2*h);
                ny = (grid[4*(i*N*N       + (j + 1)*N + k)]     - grid[4*(i*N*N       + (j - 1)*N + k)]) / (2*h);
                nz = (grid[4*(i*N*N       + j*N       + k + 1)] - grid[4*(i*N*N       + j*N       + k - 1)]) / (2*h);
                // components of normalized surface normal
                double nnx = nx / norm(nx, ny, nz);
                double nny = ny / norm(nx, ny, nz);
                double nnz = nz / norm(nx, ny, nz);

                //    // components of velocity
                //    double w1 = u1 + S * nnx;
                //    double w2 = u2 + S * nny;
                //    double w3 = u3 + S * nnz;

                // upwind finite difference approximations for partial derivatives
                double phix, phiy, phiz;

                if (w1 > 0) {
                    phix = (grid[4*(i*N*N       + j*N       + k)] - grid[4*((i - 1)*N*N + j*N       + k)]) / h;
                } else {
                    phix = (grid[4*((i + 1)*N*N + j*N       + k)] - grid[4*(i*N*N       + j*N       + k)]) / h;
                }
                if (w2 > 0) {
                    phiy = (grid[4*(i*N*N       + j*N       + k)] - grid[4*(i*N*N       + (j - 1)*N + k)]) / h;
                } else {
                    phiy = (grid[4*(i*N*N       + (j + 1)*N + k)] - grid[4*(i*N*N       + j*N       + k)]) / h;
                }
                if (w3 > 0) {
                    phiz = (grid[4*(i*N*N       + j*N       + k)] - grid[4*(i*N*N       + (j - 1)*N + k)]) / h;
                } else {
                    phiz = (grid[4*(i*N*N       + (j + 1)*N + k)] - grid[4*(i*N*N       + j*N       + k)]) / h;
                }


                // update implicit surface function
                grid[4*(i*N*N + j*N + k)] = grid[4*(i*N*N + j*N + k)]
                                                  - dt * (w1 * phix + w2 * phiy + w3 * phiz);


                // update implicit surface equation

            }
        }

    }

}

void Fire::addForce() {
    double fx, fy, fz, g;
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            for (int k = 1; k <= N; k++) {
                // compute forces
                //gravity
                g = 9.81;
                // buoyancy
                fz = alpha * (grid[4*(i*N*N + j*N + k)] - Tair);
                // vorticity
                Vector3d wijk = vort(i, j, k);
                Vector3d n((vort(i+1, j, k).norm() - wijk.norm())/h,
                        (vort(i, j+1, k).norm() - wijk.norm())/h,
                        (vort(i, j, k+1).norm() - wijk.norm())/h);
                n.normalize();
                Vector3d fconf = eps * h * n.cross(wijk);
                fconf(2) += (fz + g);
                Vector3d f = fconf;

                // add in contribution to velocity
                velNewX[i*N*N + j*N + k] = velX[i*N*N + j*N + k] + dt*f(0);
                velNewY[i*N*N + j*N + k] = velY[i*N*N + j*N + k] + dt*f(1);
                velNewZ[i*N*N + j*N + k] = velZ[i*N*N + j*N + k] + dt*f(2);
            }
        }
    }

}

Vector3d Fire::vort(int i, int j, int k) {
    double w1, w2, w3;
    w1 =   (velCZ[i*N*N     + (j+1)*N + k] - velCZ[i*N*N     + (j-1)*N + k]
          - velCY[i*N*N     + j*N + k+1]   + velCY[i*N*N     + j*N     + k-1])/(2*h);
    w2 =   (velCX[i*N*N     + j*N + k+1]   - velCX[i*N*N     + j*N     + k-1]
          - velCZ[(i+1)*N*N + j*N + k]     + velCZ[(i-1)*N*N + j*N     + k])/(2*h);
    w3 =   (velCY[(i+1)*N*N + j*N + k]     - velCY[(i-1)*N*N + j*N     + k]
          - velCX[i*N*N     + (j+1)*N + k] + velCX[i*N*N     + (j-1)*N + k])/(2*h);
    return Vector3d(w1, w2, w3);
}

void Fire::advect(vector<double> &arr, vector<double> &arrOld) {

    // solve advection

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {

                // velX is in world space units of distance / time
                // go backwards along characteristic flow line
                // --to first order this means to backwards along the velocity field
                float newI = i - dt * velX[i*N*N + j*N + k] / h;
                float newJ = j - dt * velY[i*N*N + j*N + k] / h;
                float newK = k - dt * velZ[i*N*N + j*N + k] / h;

                int x = (int) newI;
                int y = (int) newJ;
                int z = (int) newK;

                float dx = newI - x;
                float dy = newJ - y;
                float dz = newK - z;

                arr[i*N*N + j*N + k] = triLerp(x, y, z, dx, dy, dz, arrOld);

            }
        }

    }



}

double Fire::triLerp(int x, int y, int z, double dx, double dy, double dz, vector<double> &arr) {
    return  (1-dx) * (1-dy) * (1-dz) * arr[x*N*N     + y*N     + z]   +
            (1-dx) * dy     * (1-dz) * arr[x*N*N     + (y+1)*N + z]   +
            dx     * (1-dy) * (1-dz) * arr[(x+1)*N*N + y*N     + z]   +
            dx     * dy     * (1-dz) * arr[(x+1)*N*N + (y+1)*N + z]   +
            (1-dx) * (1-dy) * dz     * arr[x*N*N     + y*N     + z+1] +
            (1-dx) * dy     * dz     * arr[x*N*N     + (y+1)*N + z+1] +
            dx     * (1-dy) * dz     * arr[(x+1)*N*N + y*N     + z+1] +
            dx     * dy     * dz     * arr[(x+1)*N*N + (y+1)*N + z+1];
}


void Fire::flow() {

}


double Fire::norm(double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}

