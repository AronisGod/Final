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
    M = 1;
    ph = 0.01; // (kg/m^3) density of the "hot products"
    pf = 0.1; // (kg/m^3) density of the "fuel vapor"
    k0 = 1; // constant for dY/dt (change in time since reaction {at a point in space} over change in time) thus unit-less
    Tignition = 678; // (K) Temperature at ignition
    Tmax = 2253;  // (K) Maximum temperature
    C = S*(pf/ph - 1); // Correction for velocity discontinuity at implicit surface

    for (int n = 0; n < N*N*N; n++) {
        gridNorm.push_back(new array<double, 3>);
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
    int m = N * N * N;
    p = VectorXd(m);
    A = SparseMatrix<double>(m, m);
    typedef Triplet<double> T;
    vector<T> list;
    for (int idx = 0; idx < m; idx++) {
        list.emplace_back(idx, idx, -6);
        list.emplace_back(idx + 1, idx, 1);
        list.emplace_back(idx + N, idx, 1);
        list.emplace_back(idx + N * N, idx, 1);
        list.emplace_back(idx, idx + 1, 1);
        list.emplace_back(idx, idx + N, 1);
        list.emplace_back(idx, idx + N * N, 1);
    }
}


void Fire::propagateFront() {
    for (int n = 0; n < N*N*N; n++) {
        double nx, ny, nz;
        // components of the gradient
        nx = (grid[n + N*N] - grid[n - N*N]);
        ny = (grid[n + N]   - grid[n - N]);
        nz = (grid[n + 1]   - grid[n - 1]);
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

        if (w1 > 0) phix = (grid[n]- grid[n - N*N]) / h;
        else phix = (grid[n + N*N] - grid[n]) / h;

        if (w2 > 0) phiy = (grid[n] - grid[n - N]) / h;
        else phiy = (grid[n + N]    - grid[n]) / h;

        if (w3 > 0) phiz = (grid[n] - grid[n - 1]) / h;
        else phiz = (grid[n + 1]    - grid[n]) / h;

        // update implicit surface function
        newGrid[n] = grid[n] - dt * (w1 * phix + w2 * phiy + w3 * phiz);
    }
}

void Fire::addForce() {
    double rho, b, g = -9.81; // density, buoyancy force, accl. d.t. gravity
    for (int n = 0; n <= N*N*N; n++) {
        // compute forces
        //buoyancy
        b = alpha * (grid[n] - Tair);
        // vorticity confinement
        Vector3d wijk = vort(n);
        double Nurm = wijk.norm();
        Vector3d no((vort(n + N*N).norm() - Nurm),
                   (vort(n + N).norm()    - Nurm),
                   (vort(n + 1).norm()    - Nurm));
        no.normalize();
        Vector3d f = epsf * h * no.cross(wijk);                       // FIX DISCONTINUITY
        // Addition of vorticity confinement & buoyancy
        f(2) += b;
        // Chooses the right density
        if (grid[n] > 0) rho = pf;
        else rho = ph;                                                   // Questionable vel Choice
        // add in contribution to velocity
        velX[n] = velX[n] + dt*f(0)/rho;
        velY[n] = velY[n] + dt*f(1)/rho;
        velZ[n] = velZ[n] + dt*(f(2)/rho + g);  // plus gravity
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
    if (newGrid[n] > 0) {if (grid[n + dn] <= 0) on = -1;} // Checks two possible broadly covering
    else if (grid[n + dn] > 0) on = 1;                    // conditions for surface jumps

    if (!on) return {0,0,0};
    array<double, 3> temp = *gridNorm[n + dn];
    for (int l = 0; l < 3; l++) temp[l] = temp[l]*on*C;
    return temp;
}


void Fire::poissonPressure() {
    int NNN = N*N*N;
    VectorXd b(NNN);
    double uGrad, c = ph*h/dt;

    // Builds b vector
    for (int n = 0; n < NNN; n++) {
        // Divergence Velocity Field
        uGrad =          //Sum of (u_n+1 - u_n) for n(={x, y, z}
            velNewX[n + N*N] - velNewX[n] +
            velNewY[n + N]   - velNewY[n] +
            velNewZ[n + 1]   - velNewZ[n];

        // Handles corrections due to the discontinuity at the surface
        for (int B = 0; B < 3; B++) {
            int dn = (int) pow((float) N, 2 - B);
            if (grid[n] > 0 && grid[n + dn] <= 0)
                uGrad += C*(*gridNorm[n + dn])[B];
            else if (grid[n] <= 0 && grid[n + dn] > 0)
                uGrad += -C*(*gridNorm[n + dn])[B];
        }

        uGrad = uGrad*c;
        if (newGrid[n] > 0) uGrad = uGrad*pf/ph;
        b(n) = uGrad;
    }
    buildA();
    ConjugateGradient<SparseMatrix<double>, Lower|Upper, IncompleteCholesky<double, Lower|Upper>> cg;
    cg.compute(A);

    p = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
}


void Fire::applyPressure() {
    double rho, divP;
    for (int n = 0; n < N*N*N; n++) {
        divP = (p[n + N*N] + p[n + N] + p[n + 1] - 3*p[n]);

        if (newGrid[n] > 0) rho = pf;
        else rho = ph;

        velX[n] = velNewX[n] - dt*divP/rho;
    }
}


void Fire::updateVCenter() {
    for (int n = 0; n < N*N*N; n++) {
           // Builds velC vectors out of the vel vectors
                velCX[n] = (velX[n - N*N] + velX[n]) / 2;
                velCY[n] = (velY[n - N]   + velY[n]) / 2;
                velCZ[n] = (velZ[n - 1]   + velZ[n]) / 2;
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



Vector3d Fire::vort(int n) {
    double w1, w2, w3;

    w1 =   (velCZ[n + N]   - velCZ[n - N]   - velCY[n + 1]   + velCY[n - 1])/(2*h);
    w2 =   (velCX[n + 1]   - velCX[n - 1]   - velCZ[n + N*N] + velCZ[n - N*N])/(2*h);
    w3 =   (velCY[n + N*N] - velCY[n - N*N] - velCX[n + N]   + velCX[n - N])/(2*h);
    return Vector3d(w1, w2, w3);
}


double Fire::norm(double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}


double Fire::triLerp(int x, int y, int z, double dx, double dy, double dz, vector<double> &arr, array<double, 8> *Cor) {
    int n = x*N*N + y*N + z;
    return    (1 - dx) * (1 - dy) * (1 - dz) * (arr[n]           + (*Cor)[0]) +

              dx       * (1 - dy) * (1 - dz) * (arr[n + N*N]     + (*Cor)[4]) +
              (1 - dx) * dy       * (1 - dz) * (arr[n + N]       + (*Cor)[2]) +
              (1 - dx) * (1 - dy) * dz       * (arr[n + 1]       + (*Cor)[1]) +

              dx       * dy       * (1 - dz) * (arr[n + N*N + N]     + (*Cor)[6]) +
              dx       * (1 - dy) * dz       * (arr[n + N*N + 1]     + (*Cor)[5]) +
              (1 - dx) * dy       * dz       * (arr[n + N + 1]       + (*Cor)[3]) +

              dx       * dy       * dz       * (arr[n + N*N + N + 1] + (*Cor)[7]);
}