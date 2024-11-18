#include "grid.h"

#include <iostream>
#include <vector>

int main() {
    std::string filename = "9x9.x"; // Replace with your file path
    std::vector<double> xCoords, yCoords;

    // Read the coordinates from the file
    readCoordinates(filename, xCoords, yCoords);

    // Output the coordinates for verification
    std::cout << "X Coordinates:" << std::endl;
    for (double x : xCoords) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    std::cout << "Y Coordinates:" << std::endl;
    for (double y : yCoords) {
        std::cout << y << " ";
    }
    std::cout << std::endl;

    return 0;
}