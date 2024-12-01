#include "grid.h"
#include "solver.h"
#include "solution.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <filesystem>

/**
 * Parameters of the simulation.
 */
int n = 9; // Replace with the grid size wanted (9x9, 17x17, etc.)
double MachNumber = 0.8; // Mach number of the flow
double AoA = 0.0; // Angle of attack in degrees

/**
 * Physical properties of the flow.
 */
double rhoInf = 1.225; // Density at infinity
double pInf = 1.01325; // Pressure at infinity
double gamma = 1.4; // Heat capacity ratio
double R = 287.0; // Specific gas constant
double TInf = 288.15; // Temperature at infinity
double fluidProperties[5] = {rhoInf, pInf, gamma, R, TInf};
double dt = 0.5;
double t = 5;

/**
 * Main function.
 */
int main() {
    std::cout << "Running simulation..." << std::endl;

    // Construct the file path
    std::string filename = "../../../../NACA0012grids/" + std::to_string(n) + "x" + std::to_string(n)+".x";
    
    // Read the coordinates from the file
    std::vector<double> xCoords, yCoords; // Vectors to store the x and y coordinates of the nodes 
    readCoordinates(filename, n, xCoords, yCoords);

    // Determine the connectivity between faces and cells
    std::vector<int> faceToCellsLeft, faceToCellsRight, cellToFaces, faceToNodes;
    int faceNumber = 0, cellNumber = 0;
    connectivity(n, xCoords, yCoords, faceToCellsLeft, faceToCellsRight, cellToFaces, faceToNodes, faceNumber, cellNumber);

    // Calculate the volume of each cell
    std::vector<double> volume;
    cellVolume(n, cellNumber, xCoords, yCoords, volume);

    // Calculate the length of each face
    std::vector<double> length;
    faceLength(n, faceNumber, xCoords, yCoords, faceToNodes, length);

    // Calculate the normal vector of each face
    std::vector<double> xNormal, yNormal;
    faceNormal(n, xCoords, yCoords, xNormal, yNormal);

    // Determine the type of each face and cell
    std::vector<int> faceType, cellType;
    faceAndCellTypes(n, faceType, cellType);

    // Initialize the flow variables
    std::cout << "Initializing..." << std::endl;
    std::vector<double> W;
    std::vector<double> R;
    Initialization(n, MachNumber, AoA, fluidProperties, faceNumber, faceToCellsLeft, faceToCellsRight, xNormal, yNormal, W, R);

    // Apply the boundary conditions
    BoundaryConditions(n, MachNumber, AoA, fluidProperties, cellType, cellToFaces, xNormal, yNormal, W, R);

    // Solve
    std::cout << "Solving..." << std::endl;
    Euler(dt, t, n, MachNumber, AoA, fluidProperties, cellType, cellToFaces, faceNumber, cellNumber, volume, faceToCellsLeft, faceToCellsRight, length, xNormal, yNormal, W, R);

    // Post-process
    std::cout << "Post-processing..." << std::endl;
    // Specify the directory where the file will be created
    std::string outputDirectory = "../../../../output/";
    std::filesystem::create_directories(outputDirectory); // Create the directory if it doesn't exist
    // Write the tecplot file for rho, u, v, E and p
    std::vector<double> rho, u, v, VMag, E, p;
    calculateProperties(W, fluidProperties, cellNumber, rho, u, v, VMag, E, p);
    std::string tecplotFilename = outputDirectory + "rho.dat";
    std::string tecplotFilename2 = outputDirectory + "u.dat";
    std::string tecplotFilename3 = outputDirectory + "v.dat";
    std::string tecplotFilename4 = outputDirectory + "VMag.dat";
    std::string tecplotFilename5 = outputDirectory + "E.dat";
    std::string tecplotFilename6 = outputDirectory + "p.dat";
    WriteTecplotFile(tecplotFilename, n, n, xCoords, yCoords, rho);
    WriteTecplotFile(tecplotFilename2, n, n, xCoords, yCoords, u);
    WriteTecplotFile(tecplotFilename3, n, n, xCoords, yCoords, v);
    WriteTecplotFile(tecplotFilename4, n, n, xCoords, yCoords, VMag);
    WriteTecplotFile(tecplotFilename5, n, n, xCoords, yCoords, E);
    WriteTecplotFile(tecplotFilename6, n, n, xCoords, yCoords, p);
    // Calculate lift, drag and moment coefficients
    double cL, cD, cM; 
    aerodynamicCoefficients(n, MachNumber, AoA, fluidProperties, xCoords, yCoords, faceToCellsRight, faceToNodes, length, xNormal, yNormal, W, cL, cD, cM);
    // Write the pressure field around the airfoil
    std::string pFieldFilename = outputDirectory + "pressureField.txt";
    pressureField(n, fluidProperties, cellType, xCoords, yCoords, faceToCellsRight, faceToCellsLeft, cellToFaces, faceToNodes, W, pFieldFilename);

    return 0;
}