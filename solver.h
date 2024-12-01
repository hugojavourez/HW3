#pragma once

#include <iostream>
#include <vector>
#include <sstream>

// Initialization function
void Initialization(int n, double MachNumber, double AoA, double fluidProperties[5], int cellNumber, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R);

// Boundary condition functions
void BoundaryConditions(const int n, const double MachNumber, const double AoA, double fluidProperties[5], const std::vector<int>& cellType, std::vector<int>& cellToFaces, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R);

// Calculate the flux
void CalculateResidual( double fluidProperties[5], int faceNumber, std::vector<int> &faceToCellsLeft, std::vector<int> &faceToCellsRight, std::vector<double> &length,std::vector<double> &cellvolume, std::vector<double> &xNormal, std::vector<double> &yNormal, std::vector<double> &W,std::vector<double>& R);// rk4
void RK4(double dt,double t,  int n, double MachNumber, double AoA, double fluidProperties[5], int faceNumber, int cellNumber,std::vector<double> &cellvolume, std::vector<int> &faceToCellsLeft, std::vector<int> &faceToCellsRight, std::vector<double> &length, std::vector<double> &xNormal, std::vector<double> &yNormal, std::vector<double> &W,std::vector<double>& LC);
void WriteTecplotFile(const std::string& filename,
                      int NI, int NJ,
                      const std::vector<double>& X,
                      const std::vector<double>& Y,
                      const std::vector<double>& Volume);