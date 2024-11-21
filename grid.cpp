#include "grid.h"
#include "math.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

/**
 * Reads the coordinates from a geometry file and stores them in the provided vectors.
 *
 * @param filename The name of the file to read from.
 * @param n The grid size (nxn).
 * @param xCoords A vector to store the x-coordinates.
 * @param yCoords A vector to store the y-coordinates.
 */
void readCoordinates(const std::string& filename, const int n, std::vector<double>& xCoords, std::vector<double>& yCoords) {    
    
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
    int totalPoints = n * n;

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

/**
 * Calculates the volume of each cell in the grid.
 * 
 * @param n The grid size (nxn).
 * @param xCoords A vector containing the x-coordinates.
 * @param yCoords A vector containing the y-coordinates.
 * @param volume A vector to store the volume of each cell.
 */
void cellVolume(int n, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& volume) {
    
    int totalPoints = n * n;
    int totalCells = pow((n - 1),2);

    // Resize vectors to hold the coordinates
    volume.resize(totalCells);

    // Counter of cells
    int c = 0;

    for (int i = 0; i < totalPoints - n; i++) {
        // The node closing the layer does not contribute to a new cell
        if ((i + 1) % n != 0) {
            // Create the diagonal vectors
            double a1 = xCoords[i] - xCoords[i+n+1], a2 = yCoords[i] - yCoords[i+n+1];
            double b1 = xCoords[i+1] - xCoords[i+n], b2 = yCoords[i+1] - yCoords[i+n];
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

/**
 * Calculates the length of each face in the grid.
 * 
 * @param n The grid size (nxn).
 * @param xCoords A vector containing the x-coordinates.
 * @param yCoords A vector containing the y-coordinates.
 * @param length A vector to store the length of each face.
 */
void faceLength(int n, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& length) {
    int totalPoints = n * n;
    int totalFaces = (n - 1) * (2 * n - 1);

    // Resize vectors to hold the coordinates
    length.resize(totalFaces);

    // Counter of faces
    int f = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // The last point does not have another node to be linked to
            // And the node that close the layer has no new face (because it was computed when arrived on the layer)
            if (i + j - 2 != totalPoints || i-1 % n != 0) {
                // Create the face vector
                double A[2] = {xCoords[i+j] - xCoords[i+j+1], yCoords[i+j] - yCoords[i+j+1]};

                // Calculate the face length of the i-th cell
                norm(A,length[f]);

                // Going to the next face
                f += 1;
            }

            // The last layer of nodes does not have faces pointing to a farther layer of nodes
            // And the node that close the layer has no new face (because it was computed when arrived on the layer)
            if (j - 1 != n || i-1 % n != 0) {
                // Create the face vector
                double A[2] = {xCoords[i+j] - xCoords[i+j+n], yCoords[i+j] - yCoords[i+j+n]};

                // Calculate the face length of the i-th cell
                norm(A,length[f]);

                // Going to the next face
                f += 1;
            }
        }
    }
}

/**
 * Calculates the normal vector of each face in the grid.
 * 
 * @param n The grid size (nxn).
 * @param xCoords A vector containing the x-coordinates.
 * @param yCoords A vector containing the y-coordinates.
 * @param xNormal A vector to store the x-component of the normal vector of each face.
 * @param yNormal A vector to store the y-component of the normal vector of each face.
 */
void faceNormal(int n, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& xNormal, std::vector<double>& yNormal) {
    int totalPoints = n * n;
    int totalFaces = (n - 1) * (2 * n - 1);

    // Resize vectors to hold the coordinates
    xNormal.resize(totalFaces);
    yNormal.resize(totalFaces);

    // Counter of faces
    int f = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // The last point does not have another node to be linked to
            // And the node that close the layer has no new face (because it was computed when arrived on the layer)
            if (i + j - 2 != totalPoints || i-1 % n != 0) {
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
            if (j - 1 != n || i-1 % n != 0) {
                // Create the normal vector
                double x = yCoords[i+j] - yCoords[i+j+n];
                double y = - (xCoords[i+j] - xCoords[i+j+n]);

                // Normalize the vector
                xNormal[f] = x/(sqrt(pow(x,2)+pow(y,2)));
                yNormal[f] = y/(sqrt(pow(x,2)+pow(y,2)));

                // Going to the next face
                f += 1;
            }
        }

    }
}

/**
 * Determines the connectivity between faces and cells in the grid.
 * 
 * @param n The grid size (nxn).
 * @param faceToCellsLeft A vector to store the cells on the left of each face.
 * @param faceToCellsRight A vector to store the cells on the right of each face.
 * @param faceNumber The total number of faces.
 * @param cellNumber The total number of cells.
 */
void connectivity(int n, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, int faceNumber, int cellNumber) {
    int totalPoints = n * n;
    faceNumber = (n - 1) * (2 * n - 1);
    cellNumber = pow((n - 1),2);

    // Resize vector to hold the faces
    faceToCellsLeft.resize(faceNumber);
    faceToCellsRight.resize(faceNumber);

    // Number of the cell touching the first face of the last layer
    int lastFaces = n * (n - 2);

    for (int f = 0; f < faceNumber; f++) {
        // If the face studied is neither on the first nor last layer
        if (f > 2 * (n - 1) && f < faceNumber - n) {
            // If the face make the junction between two layers and is not the one that close a bloc
            if (f % 2 == 1 && f % (2 * (n - 1)) != 1) {
                faceToCellsLeft[f] = (f - 1) / 2;
                faceToCellsRight[f] = faceToCellsLeft[f] - 1;
            } else
            // If the face is on a layer and is not the one that close a bloc
            if (f % (2 * (n - 1)) != 1) {
                faceToCellsRight[f] = f/2;
                faceToCellsLeft[f] = faceToCellsRight[f] - (n - 1);
            } else
            // If the face closes a bloc
            {
                faceToCellsLeft[f] = (f - 2) / 2;
                faceToCellsRight[f] = faceToCellsLeft[f] + (n - 2);
            }
        } else 
        // If the face studied is on the first layer
        if (f < faceNumber - n) {
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