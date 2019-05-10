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
    int N; // Size of grid: grid will be 4 by N by N by N (4 because we have 4 quantities to keep track of)
    double dt; // size of time steps
    double h; // grid spacing
    double S; // Parameter controlling velocity of front propagation(Combustion/Reaction Rate)
    double Tair; // temperature of ambient environment
    double alpha; // Buoyancy force parameter
    double cT; // Cooling constant
    double epsh, epsf; // Vorticity confinement parameters for "hot products" & "fuel vapor" respectively
    double ph, pf; // density of the "hot products" & "fuel vapor" respectively
    double k; // constant for Y
    double Tignition; // Temperature at ignition
    double Tmax;  // maximum temperature
    double vMax;  // cap on velocity magnitude
    vector<double> grid; // grid with implicit surface at current time step
    vector<array<double, 3>*> gridNorm; // The normalized gradient field of the grid at next time step
    // We define φ to be positive in the region of space filled with fuel, negative elsewhere and zero at the reaction zone.
    vector<double> newGrid; // grid with implicit surface at next time step
    SparseMatrix<double> A;
    VectorXd p;
    vector<double> velNewX; // array with x-coordinate of velocities defined across faces of 'grid'
    vector<double> velNewY; // array with y-coordinate of velocities defined across faces of 'grid'
    vector<double> velNewZ; // array with z-coordinate of velocities defined across faces of 'grid'

    vector<double> velX; // array with x-coordinate of velocities defined across faces of 'grid'
    vector<double> velY; // array with y-coordinate of velocities defined across faces of 'grid'
    vector<double> velZ; // array with z-coordinate of velocities defined across faces of 'grid'


    // centered velocities
    vector<double> velCX; // array with x-coordinate of velocities defined across faces of 'grid'
    vector<double> velCY; // array with y-coordinate of velocities defined across faces of 'grid'
    vector<double> velCZ; // array with z-coordinate of velocities defined across faces of 'grid'

    // Y is a number such that 1- Y is the time since crossing the barrier
    vector<double> Y;
    vector<double> newY;

    // temperatures
    vector<double> T;

public:

    Fire();

    void updateY();

    void updateT();

    void propagateFront();

    double norm(double x, double y, double z);

    void step();

    void advect();

    double triLerp(int x, int y, int z, double dx, double dy, double dz, vector<double> &arr);

    void addForce();

    Vector3d vort(int i, int j, int k);

    void poissonPressure();

    void updateVCenter();

    void buildA();

};


#endif //FIRE_FIRE_H
