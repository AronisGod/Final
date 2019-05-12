#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

#ifndef FIRE_FIRE_H
#define FIRE_FIRE_H

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/IterativeLinearSolvers"

using namespace Eigen;

using namespace std;

const double INF = numeric_limits<double>::infinity();
////////////////////////////////////////////////////////

// Simulation Environment Parameters
const int N = 160; // # Cells per grid dimension
const double dt = 0.05; // (s) Time step
const double h = 0.05; // (m) Grid dimension length

// Physical Processes Parameters
const double M = 1;
const double S = 0.1; // (m/s) Reaction velocity (Combustion/Reaction Rate)
const double Tair = 300; // (K) Temperature of ambient environment
const double alpha = 0.15; // (m/Ks^2) Buoyancy force parameter
const double cT = 3000; // (K/s) Cooling constant
const double epsh = 60, epsf = 16; // Vorticity confinement parameters for "hot products" & "fuel vapor" respectively
const double ph = 0.01, pf = 0.1; // (kg/m^3) Density of the "hot products" & "fuel vapor" respectively
const double k0 = 1; // constant for Y
const double Tignition = 678; // (K) Temperature at ignition
const double Tmax = 2253;  // (K) Maximum temperature
const double vMax = h / 1.3 / dt;  // cap on velocity magnitude
const double C = S*(pf/ph - 1);     // Correction for velocity discontinuity at implicit surface

// Data Structures
// The grid is our level set.
// >0 represents fuel vapor, 0 is the reaction surface, & <=0 is gaseous products of combustion
// Each grid cell location, (i, j, k), is defined by the (x, y, z) coordinate at it's center
// Indexed into as [i*N*N + j*N + k]
array< array< array<double, N+2>, N+2>, N+2>* grid; // grid at current time step
array< array< array<double, N+2>, N+2>, N+2>* newGrid; // grid at next time step
array< array<double, 3>, N*N*N> gridNorm; // The normalized gradient field of the grid at current time step

// Pressure field p & Coefficient matrix A used in Conjugate Gradient solver for p
VectorXd p(N*N*N);
SparseMatrix<double> A;


// Velocity is defined at Cell faces before Center, e.g. At Cell(i,j,k) Vx is defined at (i-1/2, j, k)
vector<double> velX;
vector<double> velY;
vector<double> velZ;
// Velocity field after advection step, used in solving for p (above)
vector<double> velNewX;
vector<double> velNewY;
vector<double> velNewZ;
// centered velocities
vector<double> velCX;
vector<double> velCY;
vector<double> velCZ;

// Y is a number such that 1- Y is the time since crossing the barrier
vector<double> Y;
vector<double> newY;

// Temperature field
vector<double> T;


int main() {
    std::cout << "Hello, World!" << std::endl;
    Fire engine;

    return 0;
}