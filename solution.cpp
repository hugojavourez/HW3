#include "solution.h"
#include "math.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <filesystem>

/**
 * Calculates rho, u, v, E and p from the flow variables for all the domain.
 * 
 * @param W A vector containing the flow variables.
 * @param fluidProperties An array containing the fluid properties: [rhoInf, pInf, gamma, R, TInf].
 * @param cellNumber The total number of cells.
 * @param rho A vector to store the density.
 * @param u A vector to store the x-component of the velocity.
 * @param v A vector to store the y-component of the velocity.
 * @param VMag A vector to store the magnitude of the velocity.
 * @param E A vector to store the total energy.
 * @param p A vector to store the pressure.
 */
void calculateProperties(const std::vector<double>& W, const double fluidProperties[5], const int cellNumber, std::vector<double>& rho, std::vector<double>& u, std::vector<double>& v, std::vector<double>& VMag, std::vector<double>& E, std::vector<double>& p) {
    // Resize the vectors
    rho.resize(cellNumber);
    u.resize(cellNumber);
    v.resize(cellNumber);
    VMag.resize(cellNumber);
    E.resize(cellNumber);
    p.resize(cellNumber);
    
    // Loop over the cells
    for (int c = 0; c < cellNumber; c++) {
        // Calculate the fluid properties
        rho[c] = W[4 * c];
        u[c] = W[4 * c + 1] / W[4 * c];
        v[c] = W[4 * c + 2] / W[4 * c];
        VMag[c] = sqrt(pow(u[c],2) + pow(v[c],2));
        E[c] = W[4 * c + 3] / W[4 * c];
        p[c] = (fluidProperties[2] - 1) * (E[c] - 0.5 * rho[c] * (pow(u[c],2) + pow(v[c],2)));
    }
}

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

    // Print the results
    std::cout << "--- Aerodynamic coefficients ---" << std::endl;
    std::cout << "Lift coefficient: " << cL << std::endl;
    std::cout << "Drag coefficient: " << cD << std::endl;
    std::cout << "Moment coefficient: " << cM << std::endl;
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
    int cellIndex; // Index of the physical cell touching the boundary
    int otherFaceIndex; // Index of the other face that forms the cell
    double faceCenterX1, faceCenterX2; // X-coordinates of the center of the faces
    double cellCenterX; // X-coordinates of the center of the cell

    // Open the file
    std::ofstream file(pFieldFilename);

    // Loop on the faces that form the airfoil
    for (int f = (n - 1) * 4; f < (n - 1) * 6; f += 2) {
        // Find the cell on the right of the face (the physical one)
        cellIndex = faceToCellsRight[f];

        // Find the other face that forms the cell
        otherFaceIndex = cellToFaces[4*cellIndex + 3];

        // Calculate the x-coordinates of the center of the cell
        faceCenterX1 = (xCoords[faceToNodes[2 * f]] + xCoords[faceToNodes[2 * f + 1]]) / 2;
        faceCenterX2 = (xCoords[faceToNodes[2 * otherFaceIndex]] + xCoords[faceToNodes[2 * otherFaceIndex + 1]]) / 2;
        cellCenterX = (faceCenterX1 + faceCenterX2) / 2;

        // Calculate the pressure at the surface
        double pSurface = (fluidProperties[2] - 1) * ((W[4 * cellIndex + 3] / W[4 * cellIndex]) - 0.5 * W[4 * cellIndex] * pow(sqrt(pow(W[4 * cellIndex + 1],2) + pow(W[4 * cellIndex + 2],2)) / W[4 * cellIndex],2));

        // Write the pressure to the file
        file << cellCenterX << " " << pSurface << std::endl;
    }

    // Close the file
    file.close();

    // Say that the file was created successfully
    std::cout << "Pressure field file '" << pFieldFilename <<"' created successfully!" << std::endl;
}

/**
 * Writes the Tecplot file for the variable.
 * 
 * @param filename The name of the file to write.
 * @param NI
 * @param NJ
 * @param X The x-coordinates of the grid.
 * @param Y The y-coordinates of the grid.
 * @param variable The variable to write.
 */
void WriteTecplotFile(const std::string& filename, int NI, int NJ, const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& variable){
    
    int NCI = NI - 1;
    int NCJ = NJ - 1;

    std::cout << "NI: " << NI << std::endl;
    std::cout << "NJ: " << NJ << std::endl;

    std::ofstream outfile(filename);

    // Write header
    outfile << "TITLE = \"CFD Simulation Results\"\n";
    outfile << "VARIABLES = \"X\", \"Y\", \"Variable\" \n";
    outfile << "ZONE T=\"Flow Field\", I=" << NI << ", J=" << NJ - 4 << ", DATAPACKING=BLOCK\n";
    outfile << "VARLOCATION=([3-7]=CELLCENTERED)\n";

    // Crite the coordinates
    for (int idx = 2 * NI; idx < NI * NJ - 2 * NI; ++idx) {
        outfile << X[idx] << " ";
    }
    outfile << "\n";

    for (int idx = 2 * NI; idx < NI * NJ - 2 * NI; ++idx) {
        outfile << Y[idx] << " ";
    }
    outfile << "\n";

    // Write the cell variable
    auto writeCellVariable = [&](const std::vector<double>& var) {
        for (int idx = 2 * (NI - 1); idx < NCI * NCJ - 2 * (NI - 1); ++idx) {
            outfile << var[idx] << " ";
        }
        outfile << "\n";
    };
    writeCellVariable(variable);

    outfile.close();
}

/**
 * Writes the properties of the flow to several files.
 */
void writeProperties(const int n1, const int n2, const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<double>& W, const double fluidProperties[5], const int cellNumber, std::vector<double>& rho, std::vector<double>& u, std::vector<double>& v, std::vector<double>& VMag,std::vector<double>& E,std::vector<double>& p) {
    // Specify the directory where the file will be created
    std::string outputDirectory = "../../../../output/";
    std::filesystem::create_directories(outputDirectory); // Create the directory if it doesn't exist

    // Calculate the properties of the flow
    calculateProperties(W, fluidProperties, cellNumber, rho, u, v, VMag, E, p);

    // Write the tecplot file for rho, u, v, E and p
    std::string tecplotFilename = outputDirectory + "rho.dat";
    std::string tecplotFilename2 = outputDirectory + "u.dat";
    std::string tecplotFilename3 = outputDirectory + "v.dat";
    std::string tecplotFilename4 = outputDirectory + "VMag.dat";
    std::string tecplotFilename5 = outputDirectory + "E.dat";
    std::string tecplotFilename6 = outputDirectory + "p.dat";
    WriteTecplotFile(tecplotFilename, n1, n2, xCoords, yCoords, rho);
    WriteTecplotFile(tecplotFilename2, n1, n2, xCoords, yCoords, u);
    WriteTecplotFile(tecplotFilename3, n1, n2, xCoords, yCoords, v);
    WriteTecplotFile(tecplotFilename4, n1, n2, xCoords, yCoords, VMag);
    WriteTecplotFile(tecplotFilename5, n1, n2, xCoords, yCoords, E);
    WriteTecplotFile(tecplotFilename6, n1, n2, xCoords, yCoords, p);
}