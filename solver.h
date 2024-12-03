#pragma once

#include <iostream>
#include <vector>
#include <sstream>

// Initialization function
void Initialization(double MachNumber, double AoA, double fluidProperties[5], int cellNumber, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R);

// Boundary condition functions
void BoundaryConditions(const int n1, const int n2, const double MachNumber, const double AoA, double fluidProperties[5], const std::vector<int>& cellType, std::vector<int>& cellToFaces, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R);

// Calculate the flux
void CalculateResidual( double fluidProperties[5], int faceNumber, std::vector<int> &faceToCellsLeft, std::vector<int> &faceToCellsRight, std::vector<double> &length,std::vector<double> &cellvolume, std::vector<double> &xNormal, std::vector<double> &yNormal, std::vector<double> &W,std::vector<double>& R);
void RK4(double dt ,double t, int n1, int n2, double MachNumber, double AoA, double fluidProperties[5], std::vector<int>& celltype, std::vector<int>& cellToFaces, int faceNumber, int cellNumber, std::vector<double>& cellvolume, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& length, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R);
void Euler(double dt ,double t, int n1, int n2, double MachNumber, double AoA, double fluidProperties[5], std::vector<int>& celltype, std::vector<int>& cellToFaces, int faceNumber, int cellNumber, std::vector<double>& cellvolume, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& length, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R);
