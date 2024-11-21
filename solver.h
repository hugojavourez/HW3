#pragma once

#include <iostream>
#include <vector>
#include <sstream>

// Initialization function
void Initialization(const int n, const double MachNumber, const double AoA, double fluidProperties[5], std::vector<double>& W, std::vector<double>& Fc, std::vector<double>& R);

// Boundary condition functions
void BoundaryConditions(const int n, const std::vector<int>& faceType, const std::vector<int>& cellType); // Not finished
