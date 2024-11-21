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

void cellVolume(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& volume) {
    int totalPoints = n1 * n2;
    int totalCells = (n1 - 1) * (n2 - 1);

    // Resize vectors to hold the coordinates
    volume.resize(totalCells);

    // Counter of cells
    int c = 0;

    for (int i = 0; i < totalPoints - n1; i++) {
        // The node closing the layer does not contribute to a new cell
        if ((i + 1) % n2 != 0) {
            // Create the diagonal vectors
            double a1 = xCoords[i] - xCoords[i+n1+1], a2 = yCoords[i] - yCoords[i+n1+1];
            double b1 = xCoords[i+1] - xCoords[i+n1], b2 = yCoords[i+1] - yCoords[i+n1];
            double A[2] = {a1, a2};
            double B[2] = {b1, b2};
        
            // Compute the cross product between the two vectors
            crossProduct(A,B,volume[i]);

            // Calculate the volume of the i-th cell
            volume[c] = 0.5*std::abs(volume[i]);

            c += 1;
        }
    }
}

void faceLength(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& length) {
    int totalPoints = n1 * n2;
    int totalFaces = (n1 - 1) * (2 * n2 - 1);

    // Resize vectors to hold the coordinates
    length.resize(totalFaces);

    // Counter of faces
    int f = 0;

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            // The last point does not have another node to be linked to
            // And the node that close the layer has no new face (because it was computed when arrived on the layer)
            if (i + j - 2 != totalPoints || i-1 % n1 != 0) {
                // Create the face vector
                double A[2] = {xCoords[i+j] - xCoords[i+j+1], yCoords[i+j] - yCoords[i+j+1]};

                // Calculate the face length of the i-th cell
                norm(A,length[f]);

                // Going to the next face
                f += 1;
            }

            // The last layer of nodes does not have faces pointing to a farther layer of nodes
            // And the node that close the layer has no new face (because it was computed when arrived on the layer)
            if (j - 1 != n2 || i-1 % n1 != 0) {
                // Create the face vector
                double A[2] = {xCoords[i+j] - xCoords[i+j+n1], yCoords[i+j] - yCoords[i+j+n1]};

                // Calculate the face length of the i-th cell
                norm(A,length[f]);

                // Going to the next face
                f += 1;
            }
        }
    }
}

void faceNormal(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& xNormal, std::vector<double>& yNormal) {
    int totalPoints = n1 * n2;
    int totalFaces = (n1 - 1) * (2 * n2 - 1);

    // Resize vectors to hold the coordinates
    xNormal.resize(totalFaces);
    yNormal.resize(totalFaces);

    // Counter of faces
    int f = 0;

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            // The last point does not have another node to be linked to
            // And the node that close the layer has no new face (because it was computed when arrived on the layer)
            if (i + j - 2 != totalPoints || i-1 % n1 != 0) {
                // Create the normal vector
                double x = yCoords[i+j] - yCoords[i+j+1];
                double y = - (xCoords[i+j] - xCoords[i+j+1]);

                // Normalize the vector
                xNormal[f] = x/(sqrt(pow(x,2)+pow(y,2)));
                yNormal[f] = y/(sqrt(pow(x,2)+pow(y,2)));

                // Going to the next face
                f += 1;
            }

            // The last layer of nodes does not have faces pointing to a farther layer of nodes
            // And the node that close the layer has no new face (because it was computed when arrived on the layer)
            if (j - 1 != n2 || i-1 % n1 != 0) {
                // Create the normal vector
                double x = yCoords[i+j] - yCoords[i+j+n1];
                double y = - (xCoords[i+j] - xCoords[i+j+n1]);

                // Normalize the vector
                xNormal[f] = x/(sqrt(pow(x,2)+pow(y,2)));
                yNormal[f] = y/(sqrt(pow(x,2)+pow(y,2)));

                // Going to the next face
                f += 1;
            }
        }

    }
}

void connectivity(int n1, int n2, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, int faceNumber, int cellNumber) {
    int totalPoints = n1 * n2;
    faceNumber = (n1 - 1) * (2 * n2 - 1);
    cellNumber = (n1 - 1) * (n2 - 1);

    // Resize vector to hold the faces
    faceToCellsLeft.resize(faceNumber);
    faceToCellsRight.resize(faceNumber);

    // Number of the cell touching the first face of the last layer
    int lastFaces = n1 * (n2 - 2);

    for (int f = 0; f < faceNumber; f++) {
        // If the face studied is neither on the first nor last layer
        if (f > 2 * (n1 - 1) && f < faceNumber - n1) {
            // If the face make the junction between two layers and is not the one that close a bloc
            if (f % 2 == 1 && f % (2 * (n1 - 1)) != 1) {
                faceToCellsLeft[f] = (f - 1) / 2;
                faceToCellsRight[f] = faceToCellsLeft[f] - 1;
            } else
            // If the face is on a layer and is not the one that close a bloc
            if (f % (2 * (n1 - 1)) != 1) {
                faceToCellsRight[f] = f/2;
                faceToCellsLeft[f] = faceToCellsRight[f] - (n1 - 1);
            } else {
                faceToCellsLeft[f] = (f - 2) / 2;
                faceToCellsRight[f] = faceToCellsLeft[f] + (n1 - 2);
            }
        } else 
        // If the face studied is on the first layer
        if (f < faceNumber - n1) {
            faceToCellsLeft[f] = - 1;
            faceToCellsRight[f] = (f + 1) / 2;
        } else
        // If the face studied is on the last layer
        {
            faceToCellsLeft[f] = lastFaces;
            faceToCellsRight[f] = - 1;

            // Going to the next face on the last layer
            lastFaces += 1;
        }
    }
}