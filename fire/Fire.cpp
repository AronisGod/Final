//
// Created by Matthew Nicoletti on 2019-05-07.
//

#include "Fire.h"

Fire::Fire() {
    N = 160;
    h = 0.05; // (m) grid spacing
    vMax = 30;  // (m/s) cap on velocity magnitude
    dt = h / vMax / 2; // (s) size of time steps for front propagation. Fluids are updated with 5*dt
    S = 0.1; // (m/s) Parameter controlling velocity of front propagation(Combustion/Reaction Rate)
    Tair = 300; // (K) temperature of ambient environment
    alpha = 0.15; // (m/Ks^2)positive constant
    cT = 3000; // (K/s) Cooling constant
    epsh = 60;
    epsf = 16;
    ph = 0.01; // (kg/m^3) density of the "hot products"
    pf = 0.1; // (kg/m^3) density of the "fuel vapor"
    k = 1; // constant for dY/dt (change in time since reaction {at a point in space} over change in time) thus unit-less
    Tignition = 678; // (K) Temperature at ignition
    Tmax = 2253;  // (K) Maximum temperature
}

void Fire::buildA() {
    int m = N*N*N;
    p = VectorXd(m);
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

/*
double Fire::hGhost(int i, int j, int k, int c) {
    Vector3d n = *(gridNorm[i*N*N + j*N + k]);
    Vector3d v(velX[i*N*N + j*N + k], velY[i*N*N + j*N + k], velZ[i*N*N + j*N + k]);

    return v(c) - M * (1 / ph - 1 / pf) * n(c);
}

double Fire::fGhost(int i, int j, int k, int c) {
    Vector3d n = *(gridNorm[i*N*N + j*N + k]);
    Vector3d v(velX[i*N*N + j*N + k], velY[i*N*N + j*N + k], velZ[i*N*N + j*N + k]);

    return v(c) + M * (1 / ph - 1 / pf) * n(c);
}
*/

void Fire::buildA() {
    int m = N*N*N;
    p = VectorXd(m);
    A = SparseMatrix<double>(m,m);
    typedef Triplet<double> T;
    vector<T> list;
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; i < N - 1; j++) {
            for (int k = 1; i < N - 1; k++) {
                int idx = i*N*N + j*N + k;

                list.emplace_back(idx, idx, -6);
                list.emplace_back(idx + 1,   idx,       1);
                list.emplace_back(idx + N,   idx,       1);
                list.emplace_back(idx + N*N, idx,       1);
                list.emplace_back(idx,       idx + 1,   1);
                list.emplace_back(idx,       idx + N,   1);
                list.emplace_back(idx,       idx + N*N, 1);
            }
        }
    }

    A.setFromTriplets(list.begin(), list.end());
}



void Fire::propagateFront() {
    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            for (int k = 1; k < N; k++) {
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

                Vector3d nhat;
                nhat(0) = nnx;
                nhat(1) = nny;
                nhat(2) = nnz;

                gridNorm[i*N*N + j*N + k] = &nhat;

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
                velX[i*N*N + j*N + k] = velX[i*N*N + j*N + k] + dt*f(0);
                velY[i*N*N + j*N + k] = velY[i*N*N + j*N + k] + dt*f(1);
                velZ[i*N*N + j*N + k] = velZ[i*N*N + j*N + k] + dt*f(2);
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

array<double, 3> edge(int n, int dn) {
    double C = 0;
    if (grid[n] > 0 && newGrid[n] > 0) {     //Calculates Ghost values for velocities crossing barrier
        if (grid[n + dn] <= 0) C = (ph/pf - 1)*S;
    }
    else if (grid[n] <= 0 && newGrid[n] <= 0) {
        if (grid[n + dn] > 0) C = (pf/ph - 1)*S;
    }
    else if (grid[n] > 0 && newGrid[n] <= 0) {
        if (grid[n + dn] > 0) C = (pf/ph - 1)*S;
    }
    else {
        if (grid[n + dn] > 0) C = (ph/pf - 1)*S;
    }
    array<double, 3> temp = (*gridNorm[n + dn]);
    for (int l = 0; l < 3; l++) temp[l] = temp[l]*C;
    if(!C) temp = {1, 1, 1};
    return temp;
}


double Fire::triLerp(int x, int y, int z, double dx, double dy, double dz, vector<double> &arr, array<double, 8> *Cor) {
    int n = x*N*N + y*N + z;
    return    (1 - dx) * (1 - dy) * (1 - dz) * arr[n]           * (*Cor)[0] +

              dx       * (1 - dy) * (1 - dz) * arr[n + N*N]     * (*Cor)[4] +
              (1 - dx) * dy       * (1 - dz) * arr[n + N]       * (*Cor)[2] +
              (1 - dx) * (1 - dy) * dz       * arr[n + 1]       * (*Cor)[1] +

              dx       * dy       * (1 - dz) * arr[n + N*N + N] * (*Cor)[6] +
              dx       * (1 - dy) * dz       * arr[n + N*N + 1] * (*Cor)[5] +
              (1 - dx) * dy       * dz       * arr[n + N + 1]   * (*Cor)[3] +

              dx       * dy       * dz       * arr[n + N*N + N + 1] * (*Cor)[7];
}


void Fire::advect() {

        // solve advection

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {

                // velX is in world space units of 'distance / time'
                // we divide by "h" to get units of 'cells / time'
                // go backwards along characteristic flow line
                // --to first order this means to backwards along the velocity field
                int n = i*N*N + j*N + k;
                float newI = i - dt * velX[n] / h;
                float newJ = j - dt * velY[n] / h;
                float newK = k - dt * velZ[n] / h;

                    int x = (int) newI;
                    int y = (int) newJ;
                    int z = (int) newK;

                    float dx = newI - x;
                    float dy = newJ - y;
                    float dz = newK - z;

                array<array<double, 8>, 3> corrections;
                for (array<double, 8> &arr : corrections) arr.fill(0);
                array<double, 3> temp;
                for (int o = 0; o < 2; o++) {
                    for (int p = 0; p < 2; p++) {
                        for(int q = 0; q < 2; q++) {
                            int dn = o*N*N + p*N + q;
                            temp = edge(n, dn);
                            int l = o*4 + p*2 + q;
                            corrections[0][l] = temp[0];
                            corrections[1][l] = temp[1];
                            corrections[2][l] = temp[2];
                        }
                    }
                }

                velNewX[i*N*N + j*N + k] = triLerp(x, y, z, dx, dy, dz, velX, &corrections[0]);
                velNewY[i*N*N + j*N + k] = triLerp(x, y, z, dx, dy, dz, velY, &corrections[1]);
                velNewZ[i*N*N + j*N + k] = triLerp(x, y, z, dx, dy, dz, velZ, &corrections[2]);
                newY[i*N*N + j*N + k] = triLerp(x, y, z, dx, dy, dz, Y) - dt*k;
            }
        }
    }
}

    }

double Fire::norm(double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}

void Fire::poissonPressure() {
    int m = N*N*N;
    VectorXd b(m);

    double uGrad, C = ph*h/dt;
    int n;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                // Builds b vector
                n = i*N*N + j*N + k;
                uGrad =          //Sum of (u_n+1 - u_n) for n(={x, y, z}
                    velNewX[n + N*N] - velNewX[n] +
                    velNewY[n + N]   - velNewY[n] +
                    velNewZ[n + 1]   - velNewZ[n];

                // Handles corrections due to the discontinuity at the surface
                array<double, 3> temp;
                for (int B = 0; B < 3; B++) {
                    temp = edge(n, (int)pow((float)N, 2-B));
                    if(temp[0] != 1)
                        uGrad += temp[B];
                }
                // Handles particularly the case where the surface is updated across the point
                temp = edge(n, 0);
                if (temp[0] != 1) uGrad += temp[0] + temp[1] + temp[2];

                uGrad = uGrad*C;
                if (grid[n] > 0) uGrad = uGrad*pf/ph;
                b(n) = uGrad;
            }
        }
    }
    if (false) {
        buildA();
        ConjugateGradient<SparseMatrix<double>, Lower|Upper, IncompleteCholesky<double, Lower|Upper>> cg;
        cg.compute(A);
    }
    p = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
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
                newY[i * N * N + j * N + k] -= dt * k0;
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


    advect();

    poissonPressure();

    velX = velNewX;
    velY = velNewY;
    velZ = velNewZ;



    updateVCenter();

    updateT();





}


