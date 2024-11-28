#include "solution.h"
#include "math.h"

#include <iostream>
#include <vector>
#include <sstream>

void aerodynamicCoefficients(const int n, const double MachNumber, const double AoA, const double fluidProperties[5], const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<int>& faceToCellsRight, const std::vector<int>& faceToNodes, const std::vector<double>& length, const std::vector<double>& xNormal, const std::vector<double>& yNormal, const std::vector<double>& W, double& cL, double& cD, double& cM) {
    int cellIndex; // Index of the physical cell touching the boundary
    int nodeIndex1, nodeIndex2; // Index of the nodes forming the face
    double dragForce = 0, liftForce = 0; // Drag and lift forces acting on the airfoil
    double forceX, forceY, moment; // Infinitesimal forces and moments acting on the airfoil (the x and y direction are relative to the freestream velocity)
    double pSurface; // Pressure at the surface of the airfoil

    // Calculate the chord length
    double chordVector[2] = {xCoords[2 * n + int(n / 2)] - xCoords[2 * n], yCoords[2 * n + int(n / 2)] - yCoords[2 * n]}; // Vector from the first to the last node of the airfoil
    double chordLength; // Length of the chord
    norm(chordVector, chordLength);

    // Calculate the rotation matrix
    double rotationMatrix[2][2] = {{cos(AoA), -sin(AoA)}, {sin(AoA), cos(AoA)}};

    // Loop over the faces that form the airfoil
    for (int f = (n - 1) * 4; f < (n - 1) * 6; f += 2) {
        // Find the cell on the right of the face (the physical one)
        cellIndex = faceToCellsRight[f];

        // Calculate the pressure at the surface
        pSurface = (fluidProperties[2] - 1) * ((W[3 * cellIndex + 3] / W[3 * cellIndex]) - 0.5 * W[3 * cellIndex] * pow(sqrt(pow(W[3 * cellIndex + 1],2) + pow(W[3 * cellIndex + 2],2)) / W[3 * cellIndex],2));

        // Calculate the force acting on the airfoil
        forceX += pSurface * (xNormal[f] * rotationMatrix[0][0] + yNormal[f] * rotationMatrix[0][1]) * length[f];
        forceY += pSurface * (xNormal[f] * rotationMatrix[1][0] + yNormal[f] * rotationMatrix[1][1]) * length[f];
        dragForce += forceX;
        liftForce += forceY;

        // Calculate the moment acting on the airfoil
        nodeIndex1 = faceToNodes[2 * f];
        nodeIndex2 = faceToNodes[2 * f + 1];
        moment += ((xCoords[nodeIndex2] - xCoords[nodeIndex1]) * forceY - (yCoords[nodeIndex2] - yCoords[nodeIndex1]) * forceX) * length[f];
    }

    // Calculate the coefficients
    cL = liftForce / (fluidProperties[0] * pow(MachNumber,2) * chordLength);
    cD = dragForce / (fluidProperties[0] * pow(MachNumber,2) * chordLength);
    cM = moment / (fluidProperties[0] * pow(MachNumber,2) * chordLength);
}