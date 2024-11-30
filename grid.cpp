#include "grid.h"
#include "math.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <windows.h>

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
 * Determines the connectivity between nodes, faces and cells in the grid.
 * 
 * @param n The grid size (nxn).
 * @param xCoords A vector containing the x-coordinates of the nodes.
 * @param yCoords A vector containing the y-coordinates of the nodes.
 * @param faceToCellsLeft A vector to store the cells on the left of each face.
 * @param faceToCellsRight A vector to store the cells on the right of each face.
 * @param cellToFaces A vector to store the faces of each cell.
 * @param faceNumber The total number of faces.
 * @param cellNumber The total number of cells.
 */
void connectivity(const int n, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<int>& cellToFaces, std::vector<int>& faceToNodes, int& faceNumber, int& cellNumber) {
    int totalPoints = n * n;
    faceNumber = (n - 1) * (2 * n - 1);
    cellNumber = static_cast<int>(pow((n - 1), 2));

    // Resize vector to hold the faces
    faceToCellsLeft.resize(faceNumber,-2);
    faceToCellsRight.resize(faceNumber,-2);
    cellToFaces.resize(4*cellNumber,-2);
    faceToNodes.resize(2*faceNumber,-2);

    int cellIndex;
    // Counter of faces for each cell
    std::vector<int> counterFaces(cellNumber,0);

    // Index of the cell touching the first face of the last layer
    int cellIndexForLastFaces = (n - 1) * (n - 2);
    // Index of the node starting the last layer
    int nodeIndexForLastFaces = n * (n - 1);

    for (int f = 0; f < faceNumber; f++) {
        // If the face studied is neither on the first nor last layer
        if ((f >= 2 * (n - 1) - 1 || f % 2 != 0) && f < faceNumber - n + 1) {
            // If the face make the junction between two layers and is not the one that close a bloc
            if (f % 2 != 0 && f % (2 * (n - 1)) != 1) {
                faceToCellsLeft[f] = (f - 1) / 2;
                faceToCellsRight[f] = faceToCellsLeft[f] - 1;
                faceToNodes[2*f] = (f / 2) + int(f / (2 * (n - 1)));
                faceToNodes[2*f+1] = faceToNodes[2*f] + n;
            } else
            // If the face is on a layer and is not the one that close a bloc
            if (f % (2 * (n - 1)) != 1) {
                faceToCellsRight[f] = f / 2;
                faceToCellsLeft[f] = faceToCellsRight[f] - (n - 1);
                faceToNodes[2*f] = (f / 2) + int(f / (2 * (n - 1)));
                faceToNodes[2*f+1] = faceToNodes[2*f] + 1;
            } else
            // If the face closes a bloc
            {
                faceToCellsLeft[f] = f / 2;
                faceToCellsRight[f] = faceToCellsLeft[f] + (n - 2);
                faceToNodes[2*f] = n * int(f / (2 * (n - 1)));
                faceToNodes[2*f+1] = faceToNodes[2*f] + n;
            }

            // Add the face to the cells
            cellIndex = faceToCellsLeft[f];
            cellToFaces[cellIndex*4+counterFaces[cellIndex]] = f;
            counterFaces[cellIndex] += 1;

            cellIndex = faceToCellsRight[f];
            cellToFaces[cellIndex*4+counterFaces[cellIndex]] = f;
            counterFaces[cellIndex] += 1;

        } else 
        // If the face studied is on the first layer
        if (f < 2 * (n - 1)) {
            faceToCellsLeft[f] = - 1;
            faceToCellsRight[f] = (f + 1) / 2;
            faceToNodes[2*f] = f / 2;
            faceToNodes[2*f+1] = faceToNodes[2*f] + 1;

            // Add the face to the cell
            cellIndex = faceToCellsRight[f];
            cellToFaces[cellIndex*4+counterFaces[cellIndex]] = f;
            counterFaces[cellIndex] += 1;
        } else
        // If the face studied is on the last layer
        {
            faceToCellsLeft[f] = cellIndexForLastFaces;
            faceToCellsRight[f] = - 1;
            faceToNodes[2*f] = nodeIndexForLastFaces;
            faceToNodes[2*f+1] = faceToNodes[2*f] + 1;

            // Add the face to the cell
            cellIndex = faceToCellsLeft[f];
            cellToFaces[cellIndex*4+counterFaces[cellIndex]] = f;
            counterFaces[cellIndex] += 1;

            // Going to the next face and next node on the last layer
            cellIndexForLastFaces += 1;
            nodeIndexForLastFaces += 1;
        }
    }

    // Change the order of the faces for the last cell on each layer
    int faceIndex;
    for (int c = 7; c < cellNumber; c += 8) {
        faceIndex = cellToFaces[4 * c];
        cellToFaces[4 * c] = cellToFaces[4 * c + 1];
        cellToFaces[4 * c + 1] = cellToFaces[4 * c + 2];
        cellToFaces[4 * c + 2] = faceIndex;
    }
}

/**
 * Calculates the volume of each cell in the grid.
 * 
 * @param n The grid size (nxn).
 * @param cellNumber The total number of cells.
 * @param xCoords A vector containing the x-coordinates.
 * @param yCoords A vector containing the y-coordinates.
 * @param volume A vector to store the volume of each cell.
 */
void cellVolume(const int n, const int cellNumber, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& volume) {
    int totalPoints = n * n;

    // Resize vector to hold the volumes
    volume.resize(cellNumber,0.0);

    // Counter of cells
    int c = 0;

    // Loop on the nodes without the last layer
    for (int i = 0; i < totalPoints - n; i++) {
        // The node closing the layer does not contribute to a new cell
        if ((i + 1) % n != 0) {
            // Create the diagonal vectors
            double a1 = xCoords[i] - xCoords[i+n+1], a2 = yCoords[i] - yCoords[i+n+1];
            double b1 = xCoords[i+1] - xCoords[i+n], b2 = yCoords[i+1] - yCoords[i+n];
            double A[2] = {a1, a2};
            double B[2] = {b1, b2};
        
            // Compute the cross product between the two vectors
            crossProduct(A,B,volume[c]);

            // Calculate the volume of the i-th cell
            volume[c] = 0.5*std::abs(volume[c]);

            c += 1;
        }
    }
}

/**
 * Calculates the length of each face in the grid.
 * 
 * @param n The grid size (nxn).
 * @param faceNumber The total number of faces.
 * @param xCoords A vector containing the x-coordinates.
 * @param yCoords A vector containing the y-coordinates.
 * @param faceToNodes A vector containing the nodes of each face.
 * @param length A vector to store the length of each face.
 */
void faceLength(const int n, const int faceNumber, const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<int>& faceToNodes, std::vector<double>& length) {
    // Resize vector to hold the lengths
    length.resize(faceNumber,0.0);
    
    // Face vector
    double A[2];

    for (int f = 0; f < faceNumber; ++f) {
        // Get the indices of the nodes for the current face
        int node1 = faceToNodes[2 * f];
        int node2 = faceToNodes[2 * f + 1];

        // Create the face vector
        A[0] = xCoords[node1] - xCoords[node2];
        A[1] = yCoords[node1] - yCoords[node2];

        // Calculate the length of the face
        length[f] = sqrt(A[0] * A[0] + A[1] * A[1]);
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
void faceNormal(const int n, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& xNormal, std::vector<double>& yNormal) {
    int totalPoints = n * n;
    int totalFaces = (n - 1) * (2 * n - 1);

    // Resize vectors to hold the coordinates
    xNormal.resize(totalFaces,0.0);
    yNormal.resize(totalFaces,0.0);

    // Counter of faces, face vector and norm vector magnifitude
    int f = 0;
    double A[2];
    double norm;

    // Nodes index
    int current;
    int next;
    int nextLayer;

    for (int j = 0; j < n - 1; j++) { // (n-1) because the last layer of nodes does not have faces pointing to a farther layer of nodes
        for (int i = 0; i < n - 1; i++) { // (n-1) because the nodes that close the layer has no new face (because it was computed when arrived on the layer)
            // Calculate the index of the nodes
            current = i + j*n;
            next = i + j*n + 1;
            nextLayer = i + j*n + n;

            // Face on the same layer of node than the node studied
            // Create the normal vector
            A[0] = yCoords[i+j*n] - yCoords[i+j*n+1];
            A[1] = - (xCoords[i+j*n] - xCoords[i+j*n+1]);
            norm = sqrt(pow(A[0],2)+pow(A[1],2));
            // Normalize the vector
            xNormal[f] = A[0] / norm;
            yNormal[f] = A[1] / norm;
            // Going to the next face
            f += 1;

            // Face pointing to a farther layer of nodes
            // Create the normal vector
            A[0] = yCoords[i+j*n] - yCoords[i+j*n+n];
            A[1] = - (xCoords[i+j*n] - xCoords[i+j*n+n]);
            norm = sqrt(pow(A[0],2)+pow(A[1],2));
            // Normalize the vector
            xNormal[f] = A[0] / norm;
            yNormal[f] = A[1] / norm;
            // Going to the next face
            f += 1;
        }
    }
}

/**
 * Determines the type of each face and cell in the grid.
 * The number of ghost layers (for the boundary conditions) is equal to 2.
 * 
 * @param n The grid size (nxn).
 * @param faceType A vector to store the type of each face.
 * @param cellType A vector to store the type of each cell.
 * 
 * Face types:
 * -1: Wall
 *  0: Interior
 *  1: Farfield
 * 
 * Cell types:
 * -1: Wall
 *  0: Interior
 *  1: Farfield
 */
void faceAndCellTypes(const int n, std::vector<int>& faceType, std::vector<int>& cellType) {
    int totalFaces = (n - 1) * (2 * n - 1);
    int totalCells = static_cast<int>(pow((n - 1), 2));

    // Resize vectors to hold the types
    faceType.resize(totalFaces);
    cellType.resize(totalCells);

    for (int i = 0; i < totalFaces; i++) {
        if (i < 2 * 2 * (n - 1)) { // First 2 layers (Wall)
            faceType[i] = -1;
        } else if (i < totalFaces - n - 2 * (n - 1)) { // Other layers (except the last 2)
            faceType[i] = 0;
        } else { // Last 2 layers (Farfield)
            faceType[i] = 1;
        }
    }

    for (int i = 0; i < totalCells; i++) {
        if (i < 2 * (n - 1)) { // First 2 layers (Wall)
            cellType[i] = -1;
        } else if (i < totalCells - 2 * (n - 1)) { // Other layers (except the last 2)
            cellType[i] = 0;
        } else { // Last 2 layers (Farfield)
            cellType[i] = 1;
        }
    }
}