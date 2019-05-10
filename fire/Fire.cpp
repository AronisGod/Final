//
// Created by Matthew Nicoletti on 2019-05-07.
//

#include "Fire.h"





void Fire::buildA(int N) {
    int m = N*N*N;
    A = SparseMatrix<double>(m,m);
    typedef Triplet<double> T;
    vector<T> list;
    for (int idx = 0; idx < m; idx++) {
        list.emplace_back(idx, idx, -6);
        list.emplace_back(idx + 1,   idx,       1);
        list.emplace_back(idx + N,   idx,       1);
        list.emplace_back(idx + N*N, idx,       1);
        list.emplace_back(idx,       idx + 1,   1);
        list.emplace_back(idx,       idx + N,   1);
        list.emplace_back(idx,       idx + N*N, 1);
    }
    A.setFromTriplets(list.begin(), list.end());
}



void Fire::propagateFront() {
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            for (int k = 1; k <= N; k++) {
                double nx, ny, nz;
                // components of the gradient
                nx = (grid[(i+1)*N*N + j*N     + k]   - grid[(i-1)*N*N + j*N     + k])   / (2*h);
                ny = (grid[i*N*N     + (j+1)*N + k]   - grid[i*N*N     + (j-1)*N + k])   / (2*h);
                nz = (grid[i*N*N     + j*N     + k+1] - grid[i*N*N     + j*N     + k-1]) / (2*h);
                // components of normalized surface normal
                double Norm = norm(nx, ny, nz);
                double nnx = nx / Norm;
                double nny = ny / Norm;
                double nnz = nz / Norm;

                *gridNorm[i*N*N + j*N + k] = {nnx, nny, nnz};

                double w1, w2, w3;

                //    // components of velocity
                w1 = velCX[i * N * N + j * N + k] + S * nnx;
                w2 = velCY[i * N * N + j * N + k] + S * nny;
                w3 = velCZ[i * N * N + j * N + k] + S * nnz;

                // upwind finite difference approximations for partial derivatives
                double phix, phiy, phiz;

                if (w1 > 0) {
                    phix = (grid[i*N*N     + j*N     + k] - grid[(i-1)*N*N + j*N     + k]) / h;
                } else {
                    phix = (grid[(i+1)*N*N + j*N     + k] - grid[i*N*N     + j*N     + k]) / h;
                }
                if (w2 > 0) {
                    phiy = (grid[i*N*N     + j*N     + k] - grid[i*N*N     + (j-1)*N + k]) / h;
                } else {
                    phiy = (grid[i*N*N     + (j+1)*N + k] - grid[i*N*N     + j*N     + k]) / h;
                }
                if (w3 > 0) {
                    phiz = (grid[i*N*N     + j*N     + k] - grid[i*N*N     + (j-1)*N + k]) / h;
                } else {
                    phiz = (grid[i*N*N     + (j+1)*N + k] - grid[i*N*N     + j*N     + k]) / h;
                }


                // update implicit surface function
                newGrid[i*N*N + j*N + k] = grid[i*N*N + j*N + k]
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
                fz = alpha * (grid[i*N*N + j*N + k] - Tair);
                // vorticity
                Vector3d wijk = vort(i, j, k);
                Vector3d n((vort(i+1, j,   k).norm()   - wijk.norm())/h,
                           (vort(i,   j+1, k).norm()   - wijk.norm())/h,
                           (vort(i,   j,   k+1).norm() - wijk.norm())/h);
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



double Fire::norm(double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}

void Fire::poissonPressure() {
    int m = N*N*N;
    VectorXd x(m), b(m);

    double uGrad, C = ph*h/dt;
    int n;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                // Builds b vector
                n = i*N*N + j*N + k;
                uGrad =          //Sum of (u_n+1 - u_n) for n(={x, y, z}
                    velNewX[(i+1)*N*N + j*N     + k]   - velNewX[i*N*N + j*N + k] +
                    velNewY[i*N*N     + (j+1)*N + k]   - velNewY[i*N*N + j*N + k] +
                    velNewZ[i*N*N     + j*N     + k+1] - velNewZ[i*N*N + j*N + k];
                if (grid[n] > 0 && grid[n] < h*1.5) {     //Calculates Ghost values for velocities crossing barrier
                    if (grid[n + N*N] <= 0) uGrad += (ph/pf - 1)*S*(*gridNorm[n])[0];
                    if (grid[n + N] <= 0)   uGrad += (ph/pf - 1)*S*(*gridNorm[n])[1];
                    if (grid[n + 1] <= 0)   uGrad += (ph/pf - 1)*S*(*gridNorm[n])[2];
                }
                else if (grid[n] <= 0 && grid[n] > -h*1.5) {
                    if (grid[n + N*N] > 0) uGrad += (pf/ph - 1)*S*(*gridNorm[n])[0];
                    if (grid[n + N] > 0)   uGrad += (pf/ph - 1)*S*(*gridNorm[n])[1];
                    if (grid[n + 1] > 0)   uGrad += (pf/ph - 1)*S*(*gridNorm[n])[2];
                }
                uGrad = uGrad*C;
                if (grid[n] > 0) uGrad = uGrad*pf/ph;
                b(n) = uGrad;
            }
        }
    }

    ConjugateGradient<SparseMatrix<double>, Lower|Upper, IncompleteCholesky<double, Lower|Upper>> cg;
    cg.compute(A);
    x = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
// update b, and solve again
    x = cg.solve(b);
}


void Fire::updateVCenter() {

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                // Builds velC vectors out of the vel vectors
                velCX[i * N * N + j * N + k] = (velX[(i-1) * N * N + j * N + k] + velX[i * N * N + j * N + k]) / 2;
                velCY[i * N * N + j * N + k] = (velY[i * N * N + (j-1) * N + k] + velY[i * N * N + j * N + k]) / 2;
                velCZ[i * N * N + j * N + k] = (velZ[i * N * N + j * N + k-1] + velZ[i * N * N + j * N + k]) / 2;
            }
        }
    }


}

void Fire::updateY() {

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                // Builds velC vectors out of the vel vectors
                newY[i * N * N + j * N + k] -= dt * k;
            }
        }
    }
}

void Fire::updateT() {

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                // Builds velC vectors out of the vel vectors
                double Y = newY[i * N * N + j * N + k];
                if (0.9 < Y && Y < 1.0) {
                    T[i * N * N + j * N + k] = Tignition + (1 - Y) * (Tmax - Tignition) / 0.1;
                } else if (Y < 0.9) {
                    T[i * N * N + j * N + k] = Tignition;
                } else {
                    T[i * N * N + j * N + k] = Tignition + (1 - Y) * (Tmax - Tignition) / 0.1 + (1 - Y) * (1 - Y);
                }

            }
        }
    }
}



void Fire::step() {

    propagateFront();


    addForce();
    velX = velNewX;
    velY = velNewY;
    velZ = velNewZ;

    advect(velNewX, velX);
    advect(velNewY, velY);
    advect(velNewZ, velZ);


    velX = velNewX;
    velY = velNewY;
    velZ = velNewZ;

    advect(newY, Y);

    updateVCenter();

    updateT();





}


