#include "grid.h"
#include "math.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

void readCoordinates(const std::string& filename, std::vector<double>& xCoords, std::vector<double>& yCoords, int n1, int n2) {
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Read the first line (1, ignored)
    int dummy;
    file >> dummy;

    // Read the grid dimensions
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

void cellsVolume(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& volume) {
    int totalPoints = n1 * n2;

    // Resize vectors to hold the coordinates
    volume.resize(totalPoints); // totalPoints or less?

    for (int i = 0; i < totalPoints - n1; i++) {
        // Create the diagonal vectors
        double a1 = xCoords[i] - xCoords[i+n1+1], a2 = yCoords[i] - yCoords[i+n1+1];
        double b1 = xCoords[i+1] - xCoords[i+n1+2], b2 = yCoords[i+1] - yCoords[i+n1+2];
        double A[2] = {a1, a2};
        double B[2] = {b1, b2};
        
        // Compute the cross product between the two vectors
        crossProduct(A,B,volume[i]);

        // Calculate the volume of the i-th cell
        volume[i] = 0.5*std::abs(volume[i]);
    }
}

void faceLength(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& length) {
    int totalPoints = n1 * n2;

    // Resize vectors to hold the coordinates
    length.resize(totalPoints);

    for (int i = 0; i < totalPoints - n1; i++) {
        // Create the face vector
        double A[2] = {xCoords[i] - xCoords[i+1], yCoords[i] - yCoords[i+1]};

        // Calculate the face length of the i-th cell
        norm(A,length[i]);
    }
}

void showGrid()
{

}
