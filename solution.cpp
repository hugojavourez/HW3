#include "solution.h"
#include "math.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

/**
 * Calculates the aerodynamic coefficients of the airfoil.
 * 
 * @param n The grid size (nxn).
 * @param MachNumber The Mach number of the flow.
 * @param AoA The angle of attack in degrees.
 * @param fluidProperties An array containing the fluid properties: [rhoInf, pInf, gamma, R, TInf].
 * @param xCoords A vector containing the x-coordinates of the nodes.
 * @param yCoords A vector containing the y-coordinates of the nodes.
 * @param faceToCellsRight A vector containing the index of the physical cell touching the boundary.
 * @param faceToNodes A vector containing the index of the nodes forming the face.
 * @param length A vector containing the length of each face.
 * @param xNormal A vector containing the x-component of the normal vector of each face.
 * @param yNormal A vector containing the y-component of the normal vector of each face.
 * @param W A vector containing the flow variables.
 * @param cL The lift coefficient.
 * @param cD The drag coefficient.
 * @param cM The moment coefficient.
 */
void aerodynamicCoefficients(const int n, const double MachNumber, const double AoA, const double fluidProperties[5], const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<int>& faceToCellsRight, const std::vector<int>& faceToNodes, const std::vector<double>& length, const std::vector<double>& xNormal, const std::vector<double>& yNormal, const std::vector<double>& W, double& cL, double& cD, double& cM) {
    int cellIndex; // Index of the physical cell touching the boundary
    int nodeIndex1, nodeIndex2; // Index of the nodes forming the face
    double dragForce = 0, liftForce = 0; // Drag and lift forces acting on the airfoil
    double forceX = 0, forceY = 0, moment = 0; // Infinitesimal forces and moments acting on the airfoil (the x and y direction are relative to the freestream velocity)
    double pSurface; // Pressure at the surface of the airfoil

    // Calculate the chord length
    double chordVector[2] = {xCoords[2 * n + int(n / 2)] - xCoords[2 * n], yCoords[2 * n + int(n / 2)] - yCoords[2 * n]}; // Vector from the first to the last node of the airfoil
    double chordLength = 0; // Length of the chord
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
        forceX = pSurface * (xNormal[f] * rotationMatrix[0][0] + yNormal[f] * rotationMatrix[0][1]) * length[f];
        forceY = pSurface * (xNormal[f] * rotationMatrix[1][0] + yNormal[f] * rotationMatrix[1][1]) * length[f];
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

/**
 * Writes the pressure field around the airfoil in function of the x coordinates to a file.
 * 
 * @param n The grid size (nxn).
 * @param fluidProperties An array containing the fluid properties: [rhoInf, pInf, gamma, R, TInf].
 * @param cellType A vector containing the type of each cell.
 * @param xCoords A vector containing the x-coordinates of the nodes.
 * @param yCoords A vector containing the y-coordinates of the nodes.
 * @param faceToCellsRight A vector containing the index of the physical cell touching the boundary.
 * @param faceToCellsLeft A vector containing the index of the ghost cell touching the boundary.
 * @param cellToFaces A vector containing the faces of each cell.
 * @param faceToNodes A vector containing the index of the nodes forming the face.
 * @param W A vector containing the flow variables.
 * @param pFieldFilename The name of the file to write the pressure field.
 * 
 * The output file will contain the x-coordinate and the pressure at the surface of the airfoil as follows:
 * x1 p1
 * x2 p2
 * ...
 */
void pressureField(const int n, const double fluidProperties[5], const std::vector<int>& cellType, const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<int>& faceToCellsRight, const std::vector<int>& faceToCellsLeft, std::vector<int>& cellToFaces, const std::vector<int>& faceToNodes, const std::vector<double>& W, const std::string& pFieldFilename) {
    int otherFaceIndex; // Index of the other face that forms the cell
    double faceCenterX1, faceCenterX2; // X-coordinates of the center of the faces
    double cellCenterX; // X-coordinates of the center of the cell

    // Loop on the faces that form the airfoil
    for (int f = (n - 1) * 4; f < (n - 1) * 6; f += 2) {
        // Find the cell on the right of the face (the physical one)
        int cellIndex = faceToCellsRight[f];

        // Find the other face that forms the cell
        otherFaceIndex = cellToFaces[4*cellIndex + 3];

        // Calculate the x-coordinates of the center of the cell
        faceCenterX1 = (xCoords[faceToNodes[2 * f]] + xCoords[faceToNodes[2 * f + 1]]) / 2;
        faceCenterX2 = (xCoords[faceToNodes[2 * otherFaceIndex]] + xCoords[faceToNodes[2 * otherFaceIndex + 1]]) / 2;
        cellCenterX = (faceCenterX1 + faceCenterX2) / 2;

        // Calculate the pressure at the surface
        double pSurface = (fluidProperties[2] - 1) * ((W[4 * cellIndex + 3] / W[4 * cellIndex]) - 0.5 * W[4 * cellIndex] * pow(sqrt(pow(W[4 * cellIndex + 1],2) + pow(W[4 * cellIndex + 2],2)) / W[4 * cellIndex],2));

        // Write the pressure to the file
        std::ofstream file(pFieldFilename, std::ios::app);
        file << cellCenterX << " " << pSurface << std::endl;
        file.close();
    }
}