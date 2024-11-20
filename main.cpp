#include "grid.h"
#include "solver.h"
#include "math.h"

#include <iostream>
#include <vector>

int main() {
    std::string filename = "9x9.x"; // Replace with the grid size wanted
    std::vector<double> xCoords, yCoords;
    int n1, n2;

    // Read the coordinates from the file
    readCoordinates(filename, xCoords, yCoords, n1, n2);

    return 0;
}