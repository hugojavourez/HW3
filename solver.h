#pragma once

#include <iostream>
#include <vector>
#include <sstream>

// Initialization function
void Initialization(const int n, const double MachNumber, const double AoA, const double fluidProperties[5], const int faceNumber, const int cellNumber, const std::vector<int>& faceToCellsLeft, const std::vector<int>& faceToCellsRight, std::vector<double>& W, std::vector<double>& Fc, std::vector<double>& R);

// Boundary condition functions
void BoundaryConditions(const int n, const std::vector<int>& faceType, const std::vector<int>& cellType); // Not finished
