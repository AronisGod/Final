//
// Created by Matthew Nicoletti on 2019-05-07.
//

#include "Fire.h"

Fire::Fire() {
    N = 160;
    NNN =  N*N*N;
    h = 0.05; // (m) grid spacing
    vMax = 30;  // (m/s) cap on velocity magnitude
    dt = h / vMax / 2; // (s) size of time steps for front propagation. Fluids are updated with 5*dt
    S = 0.1; // (m/s) Parameter controlling velocity of front propagation(Combustion/Reaction Rate)
    Tair = 300; // (K) temperature of ambient environment
    alpha = 0.15; // (m/Ks^2)positive constant
    cT = 3000; // (K/s) Cooling constant
    epsh = 60;
    epsf = 16;
    M = 1;
    ph = 0.01; // (kg/m^3) density of the "hot products"
    pf = 0.1; // (kg/m^3) density of the "fuel vapor"
    k0 = 1; // constant for dY/dt (change in time since reaction {at a point in space} over change in time) thus unit-less
    Tignition = 678; // (K) Temperature at ignition
    Tmax = 2253;  // (K) Maximum temperature
    C = S*(pf/ph - 1); // Correction for velocity discontinuity at implicit surface

    for (int n = 0; n < NNN; n++)
        gridNorm.push_back(new array<double, 3>);

    for (int n = 0; n < (int)pow(N+1, 3); n++){
        grid.push_back(INF);
        newGrid.push_back(INF);

    }
}


void Fire::step() {

    propagateFront();

    addForce();

    advect();

    poissonPressure();

    applyPressure();

    updateVCenter();

    updateT();
}


void Fire::buildA() {
    p = VectorXd(NNN);
    A = SparseMatrix<double>(NNN, NNN);
    typedef Triplet<double> T;
    vector<T> list;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int idx = i*N*N + j*N + k;

                list.emplace_back(idx, idx, -6);
                if (k != N-1) list.emplace_back(idx, idx + 1, 1);
                if (j != N-1) list.emplace_back(idx, idx + N, 1);
                if (i != N-1) list.emplace_back(idx, idx + N * N, 1);

                if (k != 0) list.emplace_back(idx, idx - 1, 1);
                if (j != 0) list.emplace_back(idx, idx - N, 1);
                if (i != 0) list.emplace_back(idx, idx - N * N, 1);
            }
        }
    }
}


void Fire::propagateFront() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int n = i*N*N + j*N + k*N;
                double nx, ny, nz;
                // components of the gradient
                if      (i == 0)   nx = (grid[n + N*N] - grid[n])*2;
                else if (i == N-1) nx = (grid[n] - grid[n - N*N])*2;
                else               nx = (grid[n + N*N] - grid[n - N*N]);
                ////////////////////////////
                if      (j == 0)   ny = (grid[n + N] - grid[n])*2;
                else if (j == N-1) ny = (grid[n] - grid[n - N])*2;
                else               ny = (grid[n + N] - grid[n - N]);
                ////////////////////////////
                if      (k == 0)   nz = (grid[n + 1] - grid[n])*2;
                else if (k == N-1) nz = (grid[n] - grid[n - 1])*2;
                else               nz = (grid[n + 1] - grid[n - 1]);
                // components of normalized surface normal
                double Norm = norm(nx, ny, nz);
                nx = nx / Norm;
                ny = ny / Norm;
                nz = nz / Norm;

                *gridNorm[n] = {nx, ny, nz};

                double w1, w2, w3;

                //    // components of velocity
                w1 = velCX[n] + S * nx;
                w2 = velCY[n] + S * ny;
                w3 = velCZ[n] + S * nz;

                // upwind finite difference approximations for partial derivatives
                double phix, phiy, phiz;

                if (i == N-1 || (i != 0 && w1 > 0)) phix = (grid[n] - grid[n - N * N]) / h;
                else phix = (grid[n + N * N] - grid[n]) / h;

                if (j == N-1 || (j != 0 && w2 > 0)) phiy = (grid[n] - grid[n - N]) / h;
                else phiy = (grid[n + N] - grid[n]) / h;

                if (k == N-1 || (k != 0 && w3 > 0)) phiz = (grid[n] - grid[n - 1]) / h;
                else phiz = (grid[n + 1] - grid[n]) / h;

                // update implicit surface function
                newGrid[n] = grid[n] - dt * (w1 * phix + w2 * phiy + w3 * phiz);
            }
        }
    }
}

void Fire::addForce() {
    double rho, b, g = -9.81; // density, buoyancy force, accl. d.t. gravity
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int n = i*N*N + j*N + k;
                // compute forces
                //buoyancy
                b = alpha * (grid[n] - Tair);
                // vorticity confinement
                Vector3d wijk = vort(i, j, k);
                double Nurm = wijk.norm();
                Vector3d no((vort(i+1, j, k).norm() - Nurm),
                            (vort(i, j+1, k).norm() - Nurm),
                            (vort(i, j, k+1).norm() - Nurm));
                no.normalize();
                Vector3d f = epsf * h * no.cross(wijk);                       // FIX DISCONTINUITY
                // Addition of vorticity confinement & buoyancy
                f(2) += b;
                // Chooses the right density
                if (grid[n] > 0) rho = pf;
                else rho = ph;                                                   // Questionable vel Choice
                // add in contribution to velocity
                velX[n] = velX[n] + dt * f(0) / rho;
                velY[n] = velY[n] + dt * f(1) / rho;
                velZ[n] = velZ[n] + dt * (f(2) / rho + g);  // plus gravity
            }
        }
    }
}


void Fire::advect() {
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

                    int x = clamp((int) newI);
                    int y = clamp((int) newJ);
                    int z = clamp((int) newK);


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
                velNewX[n] = triLerp(x, y, z, dx, dy, dz, velX, &corrections[0]);
                velNewY[n] = triLerp(x, y, z, dx, dy, dz, velY, &corrections[1]);
                velNewZ[n] = triLerp(x, y, z, dx, dy, dz, velZ, &corrections[2]);
                newY[n]    = triLerp(x, y, z, dx, dy, dz, Y) - dt*k;
            }
        }
    }
}


array<double, 3> Fire::edge(int n, int dn) {
    int on = 0;
    if (n + dn > N-1) on = 0;    // Boundary Check
    else if (newGrid[n] > 0) {if (grid[n + dn] <= 0) on = -1;} // Checks two possible broadly covering
    else if (grid[n + dn] > 0) on = 1;                    // conditions for surface jumps

    if (!on) return {0,0,0};
    array<double, 3> temp = *gridNorm[n + dn];
    for (int l = 0; l < 3; l++) temp[l] = temp[l]*on*C;
    return temp;
}


void Fire::poissonPressure() {
    VectorXd b(NNN);
    double uDiv, c = ph*h/dt;

    // Builds b vector
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int n = i*N*N + j*N + k;
                uDiv = 0;
                // Divergence Velocity Field
                //Sum of (u_n+1 - u_n) for n(={x, y, z}
                if (i == N-1) uDiv += -velNewX[n];
                else    uDiv += velNewX[n + N*N] - velNewX[n];
                if (i == N-1) uDiv += -velNewX[n];
                else    uDiv += velNewY[n + N] - velNewY[n];
                if (i == N-1) uDiv += -velNewX[n];
                else    uDiv += velNewZ[n + 1] - velNewZ[n];

                // Handles corrections due to the discontinuity at the surface
                for (int B = 0; B < 3; B++) {
                    int dn = (int) pow((float) N, 2 - B);
                    if (n + dn > N-1) continue;
                    if (grid[n] > 0 && grid[n + dn] <= 0)
                        uDiv += C * (*gridNorm[n + dn])[B];
                    else if (grid[n] <= 0 && grid[n + dn] > 0)
                        uDiv += -C * (*gridNorm[n + dn])[B];
                }

                uDiv = uDiv * c;
                if (newGrid[n] > 0) uDiv = uDiv * pf / ph;
                b(n) = uDiv;
            }
        }
    }
    buildA();
    ConjugateGradient<SparseMatrix<double>, Lower|Upper, IncompleteCholesky<double, Lower|Upper>> cg;
    cg.compute(A);

    p = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
}


void Fire::applyPressure() {
    double rho, vMag;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int n = i*N*N + j*N + k;
                if (newGrid[n] > 0) rho = pf;
                else rho = ph;

                if (i == 0) velX[n] = velNewX[n] - dt * (p[n + N*N] - p[n])*2 / rho;
                else if (i == N-1) velX[n] = velNewX[n] - dt * (p[n] - p[n - N*N])*2 / rho;
                else velX[n] = velNewX[n] - dt * (p[n + N*N] - p[n - N*N]) / rho;

                if (j == 0) velY[n] = velNewY[n] - dt * (p[n + N] - p[n])*2 / rho;
                else if (j == N-1) velY[n] = velNewY[n] - dt * (p[n] - p[n - N])*2 / rho;
                else velY[n] = velNewY[n] - dt * (p[n + N] - p[n - N]) / rho;

                if (i == 0) velZ[n] = velNewZ[n] - dt * (p[n + 1] - p[n])*2 / rho;
                else if (i == N-1) velZ[n] = velNewZ[n] - dt * (p[n] - p[n - 1])*2 / rho;
                else velZ[n] = velNewZ[n] - dt * (p[n + 1] - p[n - 1]) / rho;


                vMag = norm(velX[n], velY[n], velZ[n]);
                // Clamps the final velocity to a maximum
                if (vMag > vMax) {
                    velX[n] = vMax * velX[n] / vMag;
                    velY[n] = vMax * velY[n] / vMag;
                    velZ[n] = vMax * velZ[n] / vMag;
                }
            }
        }
    }
}


void Fire::updateVCenter() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                int n = i*N*N + j*N + k;
                // Builds velC vectors out of the vel vectors
                if (i == 0) velCX[n] = velX[n] / 2;
                else velCX[n] = (velX[n - N*N] + velX[n]) / 2;
                if (j == 0) velCY[n] = velY[n] / 2;
                else velCY[n] = (velY[n - N] + velY[n]) / 2;
                if (k == 0) velCZ[n] = velZ[n] / 2;
                else velCZ[n] = (velZ[n - 1] + velZ[n]) / 2;
            }
        }
    }
}


void Fire::updateT() {
    for (int n = 0; n < N*N*N; n++) {
        // Builds velC vectors out of the vel vectors
        double Y = newY[n];
        if (0.9 < Y && Y < 1.0)
            T[n] = Tignition + (1 - Y) * (Tmax - Tignition) / 0.1;
        else if (Y < 0.9)
            T[n] = Tignition;
        else
            T[n] = Tignition + (1 - Y) * (Tmax - Tignition) / 0.1 + (1 - Y) * (1 - Y);
    }
}



Vector3d Fire::vort(int i, int j, int k) {
    double w1, w2, w3;
    int n = i*N*N + j*N + k;
    if (i < 0 || i == N || j < 0 || j == N || k < 0 || k == N)
        return Vector3d(0, 0, 0);
    double Xy, Xz, Yx, Yz, Zx, Zy;
    if (i == 0)       {Yx =  velCY[n + N*N]; Zx =  velCZ[n + N*N];}
    else if(i == N-1) {Yx = -velCY[n - N*N]; Zx = -velCZ[n - N*N];}
    else {Yx = velCY[n + N*N] - velCY[n - N*N]; Zx = velCZ[n + N*N] - velCZ[n - N*N];}

    if (j == 0)       {Xy =  velCX[n + N]; Zy =  velCZ[n + N];}
    else if(j == N-1) {Xy = -velCX[n - N]; Zy = -velCZ[n - N];}
    else {Xy = velCX[n + N] - velCX[n - N]; Zy = velCZ[n + N] - velCZ[n - N];}

    if (k == 0)       {Xz =  velCX[n + 1]; Yz =  velCY[n + 1];}
    else if(k == N-1) {Xz = -velCX[n - 1]; Yz = -velCY[n - 1];}
    else {Xz = velCX[n + 1] - velCX[n - 1]; Yz = velCY[n + 1] - velCY[n - 1];}

    w1 =   (Zy - Yz)/(2*h);
    w2 =   (Xz - Zx)/(2*h);
    w3 =   (Yx - Xy)/(2*h);
    return Vector3d(w1, w2, w3);
}


double Fire::norm(double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}


double Fire::triLerp(int x, int y, int z, double dx, double dy, double dz, vector<double> &arr, array<double, 8> *Cor) {
    int n = x*N*N + y*N + z; // For Convenience
    // Boundary Checks
    if (x == N-1) dx = 0;
    else if (x == 0) dx = 1;
    if (y == N-1) dy = 0;
    else if (y == 0) dy = 1;
    if (z == N-1) dz = 0;
    else if (z == 0) dz = 1;

    return    (1 - dx) * (1 - dy) * (1 - dz) * (arr[n]           + (*Cor)[0]) +

              dx       * (1 - dy) * (1 - dz) * (arr[n + N*N]     + (*Cor)[4]) +
              (1 - dx) * dy       * (1 - dz) * (arr[n + N]       + (*Cor)[2]) +
              (1 - dx) * (1 - dy) * dz       * (arr[n + 1]       + (*Cor)[1]) +

              dx       * dy       * (1 - dz) * (arr[n + N*N + N]     + (*Cor)[6]) +
              dx       * (1 - dy) * dz       * (arr[n + N*N + 1]     + (*Cor)[5]) +
              (1 - dx) * dy       * dz       * (arr[n + N + 1]       + (*Cor)[3]) +

              dx       * dy       * dz       * (arr[n + N*N + N + 1] + (*Cor)[7]);
}