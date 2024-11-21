#include "grid.h"
#include "solver.h"
#include "math.h"

#include <iostream>
#include <vector>

int n = 9; // Replace with the grid size wanted (9x9, 17x17, etc.)

int main() {

    std::string filename = std::to_string(n) + "x" + std::to_string(n) + ".txt"; // Name of the geometry file to read from
    std::vector<double> xCoords, yCoords; // Vectors to store the x and y coordinates of the nodes

    // Read the coordinates from the file
    readCoordinates(filename, n, xCoords, yCoords);

    return 0;
}