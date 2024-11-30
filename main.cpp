#include "grid.h"
#include "solver.h"
#include "solution.h"
#include "test.h"

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
double pInf = 101325.0; // Pressure at infinity
double gamma = 1.4; // Heat capacity ratio
double R = 287.0; // Specific gas constant
double TInf = 288.15; // Temperature at infinity
double fluidProperties[5] = {rhoInf, pInf, gamma, R, TInf};

/**
 * Main function.
 */
int main() {
    // Get the file path of the grid
    std::string filename = "../../../../NACA0012grids/" + std::to_string(n) + "x" + std::to_string(n) + ".x"; // Construct the file path

    // Read the coordinates from the file
    std::vector<double> xCoords, yCoords;
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

    // Specify the directory where the file will be created
    std::string outputDirectory = "../../../../output/";
    std::filesystem::create_directories(outputDirectory); // Create the directory if it doesn't exist

    // Tests
    std::cout << "Running tests..." << std::endl;
    std::string tecplotFilename = outputDirectory + "volumeTest.dat"; // Name of the file to write the tecplot file
    volumeTest(n, cellNumber, xCoords, yCoords, cellToFaces, faceToNodes, volume, tecplotFilename);


    // Initialize the flow variables
    std::vector<double> W;
    std::vector<double> R;
    Initialization(n, MachNumber,  AoA,  fluidProperties,  faceNumber, faceToCellsLeft, faceToCellsRight, xNormal, yNormal, W, R);

    // Apply the boundary conditions
    BoundaryConditions(n, MachNumber,  AoA, fluidProperties, cellType, cellToFaces, xNormal, yNormal, W, R);

    // Solve
    std::vector<double> Fc;
    // ...

    // Post-process
    double cL, cD, cM; // Lift, drag and moment coefficients
    aerodynamicCoefficients(n, MachNumber, AoA, fluidProperties, xCoords, yCoords, faceToCellsRight, faceToNodes, length, xNormal, yNormal, W, cL, cD, cM);
    std::string pFieldFilename = outputDirectory + "pressureField.txt"; // Name of the file to write the pressure field around the airfoil
    pressureField(n, fluidProperties, cellType, xCoords, yCoords, faceToCellsRight, faceToCellsLeft, cellToFaces, faceToNodes, W, pFieldFilename);

    return 0;
}