#include "grid.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

void readCoordinates(const std::string& filename, std::vector<double>& xCoords, std::vector<double>& yCoords) {
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Read the first line (1, ignored)
    int dummy;
    file >> dummy;

    // Read the grid dimensions
    int n1, n2;
    file >> n1 >> n2;

    // Total number of points
    int totalPoints = n1 * n2;

    // Resize vectors to hold the coordinates
    xCoords.resize(totalPoints);
    yCoords.resize(totalPoints);

    // Read x-coordinates
    for (int i = 0; i < totalPoints; ++i) {
        file >> xCoords[i];
    }

    // Read y-coordinates
    for (int i = 0; i < totalPoints; ++i) {
        file >> yCoords[i];
    }

    file.close();
}

void showGrid()
{
    
}
