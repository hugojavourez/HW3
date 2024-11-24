#pragma once

#include <iostream>
#include <vector>
#include <sstream>

// Initialization function
void Initialization(int n, double MachNumber, double AoA, double fluidProperties[5], int faceNumber, int cellNumber, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& Fc, std::vector<double>& R);

// Boundary condition functions
void BoundaryConditions(const int n, const std::vector<int>& faceType, const std::vector<int>& cellType); // Not finished
void Solve();
