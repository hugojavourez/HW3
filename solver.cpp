#include "math.h"

#include <iostream>
#include <vector>

void Initialization(const int n, const double MachNumber, const double AoA, double fluidProperties[5], std::vector<double>& W, std::vector<double>& Fc, std::vector<double>& R) {
    std::cout << "Initialization function called" << std::endl;

    // Calculate the flow velocity in each direction
    double flowDirection[2] = {std::cos(AoA), std::sin(AoA)};
    double flowVelocity[2] = {flowDirection[0] * MachNumber, flowDirection[1] * MachNumber};


}