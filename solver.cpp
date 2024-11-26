#include "math.h"

#include <iostream>
#include <vector>

void Initialization(int n, double MachNumber, double AoA, double fluidProperties[5], int cellNumber, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R) {
    std::cout << "Initialization function called" << std::endl;

    // Rezise vectors to hold the values
    W.resize(4 * cellNumber);
    R.resize(cellNumber);

    // Calculate the flow velocity in each direction
    double flowDirection[2] = {cos(AoA), sin(AoA)};
    double flowVelocity[2] = {flowDirection[0] * MachNumber, flowDirection[1] * MachNumber};
    double velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));

    // Calculate the total energy
    double totalEnergy = (fluidProperties[1] / ((fluidProperties[2] - 1)) + 0.5 * fluidProperties[0] * pow(velocityMagnitude,2));

    // Calculate all vectors
    for (int c = 0; c < cellNumber; c++) {
        // Calculate the fluid properties vector W
        W[3 * c] = fluidProperties[0];
        W[3 * c + 1] = fluidProperties[0] * flowVelocity[0];
        W[3 * c + 2] = fluidProperties[0] * flowVelocity[1];
        W[3 * c + 3] = fluidProperties[0] * totalEnergy;

        // Calculate the residual vector R (initially set to zero)
        R[c] = 0;
    }
}

void Solve(){

}