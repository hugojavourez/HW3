#include "grid.h"
#include "solver.h"
#include "math.h"

#include <iostream>
#include <vector>

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

    std::string filename = "NACA0012grids/" + std::to_string(n) + "x" + std::to_string(n) + ".x"; // Construct the file path
    
    std::vector<double> xCoords, yCoords; // Vectors to store the x and y coordinates of the nodes

    // Read the coordinates from the file
    readCoordinates(filename, n, xCoords, yCoords);

    // Determine the connectivity between faces and cells
    std::vector<int> faceToCellsLeft, faceToCellsRight, cellToFaces, faceToNodes;
    int faceNumber = 0, cellNumber = 0;
    connectivity(n, xCoords, yCoords, faceToCellsLeft, faceToCellsRight, cellToFaces, faceToNodes, faceNumber, cellNumber);

    // Calculate the volume of each cell
    std::vector<double> volume;
    cellVolume(n, xCoords, yCoords, volume);

    // Calculate the length of each face
    std::vector<double> length;
    faceLength(n, xCoords, yCoords, length);

    // Calculate the normal vector of each face
    std::vector<double> xNormal, yNormal;
    faceNormal(n, xCoords, yCoords, xNormal, yNormal);

    // Determine the type of each face and cell
    std::vector<int> faceType, cellType;
    faceAndCellTypes(n, faceType, cellType);

    // Initialize the flow variables
    std::vector<double> W;
    std::vector<double> R;
    Initialization(n, MachNumber,  AoA,  fluidProperties,  faceNumber, faceToCellsLeft, faceToCellsRight, xNormal, yNormal, W, R);

    // Apply the boundary conditions
    BoundaryConditions(n, MachNumber,  AoA, fluidProperties, cellType, cellToFaces, xNormal, yNormal, W, R);

    // Solve
    std::vector<double> Fc;
    // ...
    
    return 0;
}