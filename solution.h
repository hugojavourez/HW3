#pragma once

#include <iostream>
#include <vector>
#include <sstream>

void calculateProperties(const std::vector<double>& W, const double fluidProperties[5], const int cellNumber, std::vector<double>& rho, std::vector<double>& u, std::vector<double>& v, std::vector<double>& VMag, std::vector<double>& E, std::vector<double>& p);
void aerodynamicCoefficients(const int n, const double MachNumber, const double AoA, const double fluidProperties[5], const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<int>& faceToCellsRight, const std::vector<int>& faceToNodes, const std::vector<double>& length, const std::vector<double>& xNormal, const std::vector<double>& yNormal, const std::vector<double>& W, double& cL, double& cD, double& cM);
void techplot(const int n, const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<int>& faceToNodes, const std::vector<double>& W, const std::string& techplotFilename);
void pressureField(const int n, const double fluidProperties[5], const std::vector<int>& cellType, const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<int>& faceToCellsRight, const std::vector<int>& faceToCellsLeft, std::vector<int>& cellToFaces, const std::vector<int>& faceToNodes, const std::vector<double>& W, const std::string& pFieldFilename);
