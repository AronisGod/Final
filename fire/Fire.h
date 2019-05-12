//
// Created by Matthew Nicoletti on 2019-05-07.
//
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>

#ifndef FIRE_FIRE_H
#define FIRE_FIRE_H

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/IterativeLinearSolvers"

using namespace Eigen;

using namespace std;


class Fire {
private:
    // Simulation Environment Parameters
    int N; // Number of cells per grid axis
    double dt; // size of time steps
    double h; // Grid cell edge length

    // Physical Processes Parameters
    double M;
    double S; // Parameter controlling velocity of front propagation(Combustion/Reaction Rate)
    double Tair; // temperature of ambient environment
    double alpha; // Buoyancy force parameter
    double cT; // Cooling constant
    double epsh, epsf; // Vorticity confinement parameters for "hot products" & "fuel vapor" respectively
    double ph, pf; // density of the "hot products" & "fuel vapor" respectively
    double k0; // constant for Y
    double Tignition; // Temperature at ignition
    double Tmax;  // maximum temperature
    double vMax;  // cap on velocity magnitude
    double C;     // Correction for velocity discontinuity at implicit surface

    // Data Structures
    // The grid is our level set.
    // >0 represents fuel vapor, 0 is the reaction surface, & <=0 is gaseous products of combustion
    // Each grid cell location, (i, j, k), is defined by the (x, y, z) coordinate at it's center
    // Indexed into as [i*N*N + j*N + k]
    vector<double> grid; // grid at current time step
    vector<double> newGrid; // grid at next time step
    vector< array<double, 3>* > gridNorm; // The normalized gradient field of the grid at current time step

    // Pressure field p & Coefficient matrix A used in Conjugate Gradient solver for p
    VectorXd p;
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

public:
    //Constructor
    Fire();

    // Simulation Terms
    void propagateFront();   // Moves the reaction front forward in time

    void addForce();         // Gravity/Buoyancy/Vorticity Force Terms

    void advect();           // Advection

    void poissonPressure();  // Calculates pressure field which upholds incompressibility

    void applyPressure();    // Applies the newly found pressure field to update our final velocity of the step

    void updateT();          // Updates Temperature field based on Y data (Reaction time history)

    void step();             // Calculates one entire frame

    // Sub calculations
    array<double, 3> edge(int n, int dn); // Returns correction to velocity at n + dn if it crosses the reaction boundary

    Vector3d vort(int n);    // Returns the curl of Velocity at grid index n AKA the vorticity

    void updateVCenter();    // Updates Velocity field defined at cell centers

    void buildA();           // Constructs the coefficient matrix used to find pressure field via Conjugate Gradient

    // Helper Functions
    double triLerp(int x, int y, int z, double dx, double dy, double dz, vector<double> &arr, array<double, 8>* = NULL); // Tri-linear Interpolation

    double norm(double x, double y, double z);        // magnitude of vector with given components

};


#endif //FIRE_FIRE_H
