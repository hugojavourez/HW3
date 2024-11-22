#include "math.h"

#include <iostream>
#include <vector>

void Initialization(const int n, const double MachNumber, const double AoA, const double fluidProperties[5], const int faceNumber, const int cellNumber, const std::vector<int>& faceToCellsLeft, const std::vector<int>& faceToCellsRight, const std::vector<double>& xNormal, const std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& Fc, std::vector<double>& R) {
    std::cout << "Initialization function called" << std::endl;

    // Rezise vectors to hold the values
    W.resize(4 * faceNumber);
    Fc.resize(4 * faceNumber);
    R.resize(cellNumber);

    // Calculate the flow velocity in each direction
    double flowDirection[2] = {cos(AoA), sin(AoA)};
    double flowVelocity[2] = {flowDirection[0] * MachNumber, flowDirection[1] * MachNumber};
    double velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));

    // Initialize the contravariant velocity and the normal vector
    double Vc;
    double normalVector[2];

    // Calculate the total energy
    double totalEnergy = (fluidProperties[1] / ((fluidProperties[2] - 1)) + 0.5 * fluidProperties[0] * pow(velocityMagnitude,2));

    // Calculate all vectors
    for (int f = 0; f < faceNumber; f++) {
        // Calculate W
        W[3 * f] = fluidProperties[0];
        W[3 * (f + 1)] = fluidProperties[0] * flowVelocity[0];
        W[3 * (f + 2)] = fluidProperties[0] * flowVelocity[1];
        W[3 * (f + 3)] = fluidProperties[0] * totalEnergy;

        // Calculate the contravariant velocity
        normalVector[0] = xNormal[f];
        normalVector[1] = yNormal[f];
        dotProduct(flowVelocity, normalVector, Vc);

        // Calculate Fc
        Fc[3 * f] = fluidProperties[0] * Vc;
        Fc[3 * (f + 1)] = fluidProperties[0] * Vc * flowVelocity[0] + fluidProperties[1] * xNormal[f];
        Fc[3 * (f + 2)] = fluidProperties[0] * Vc * flowVelocity[1] + fluidProperties[1] * yNormal[f];
        Fc[3 * (f + 3)] = fluidProperties[0] * Vc * totalEnergy + fluidProperties[1] * Vc;

        // Calculate the residual vector R
        R[faceToCellsLeft[f]] += 0;
        R[faceToCellsRight[f]] -= 0;
    }
}