#include "test.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>

/**
 * Creates a tecplot file with the cell centroids and volumes.
 * 
 * @param n The grid size (nxn).
 * @param cellNumber The number of cells.
 * @param xCoords A vector containing the x-coordinates of the nodes.
 * @param yCoords A vector containing the y-coordinates of the nodes.
 * @param cellToFaces A vector containing the faces of each cell.
 * @param faceToNodes A vector containing the nodes of each face.
 * @param volume A vector containing the volume of each cell.
 * @param filename The name of the tecplot file to create.
 */
void volumeTest(const int n, const int cellNumber, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<int>& cellToFaces, std::vector<int>& faceToNodes, std::vector<double>& volume, std::string filename) {
    double xCoordCell, yCoordCell; // Coordinates of the cell centroid
    double xCoordFace1, yCoordFace1; // Coordinates of the face 1 centroid
    double xCoordFace2, yCoordFace2; // Coordinates of the face 2 centroid

    // Open the tecplot file
    std::ofstream file(filename);

    // Write the header of the tecplot file
    file << "TITLE = \"Cell centroids and volumes\"" << std::endl;
    file << "VARIABLES = \"X\", \"Y\", \"Volume\"" << std::endl;
    file << "ZONE T=\"Volume\", N=" << n * n << ", E=" << cellNumber <<", F=FEPOINT" << std::endl;

    // Write the node coordinates
    for (int i = 0; i < n * n; i++) {
        file << std::fixed << std::setprecision(6) << xCoords[i] << " " << yCoords[i] << std::endl;
    }

    // Write the element connectivity (the nodes forming each cell)
    for (int i = 0; i < cellNumber; i++) {
        file << faceToNodes[cellToFaces[i*4]*2]+1 << " "
             << faceToNodes[cellToFaces[i*4]*2+1]+1 << " "
             << faceToNodes[cellToFaces[i*4+3]*2+1]+1 << " "
             << faceToNodes[cellToFaces[i*4+3]*2]+1 << std::endl;
    }

    // Loop on the cells
    for (int i = 0; i < cellNumber; i++) {
        // Get the coordinates of the faces (which form the cell) centroids
        // Get the coordinates of the face 1 centroid
        xCoordFace1 = 0.5*(xCoords[faceToNodes[cellToFaces[i*4]*2]] + xCoords[faceToNodes[cellToFaces[i*4]*2+1]]);
        yCoordFace1 = 0.5*(yCoords[faceToNodes[cellToFaces[i*4]*2]] + yCoords[faceToNodes[cellToFaces[i*4]*2+1]]);
        // Get the coordinates of the face 2 centroid
        xCoordFace2 = 0.5*(xCoords[faceToNodes[cellToFaces[i*4+3]*2]] + xCoords[faceToNodes[cellToFaces[i*4+3]*2+1]]);
        yCoordFace2 = 0.5*(yCoords[faceToNodes[cellToFaces[i*4+3]*2]] + yCoords[faceToNodes[cellToFaces[i*4+3]*2+1]]);

        // Get the coordinates of the cell centroid
        xCoordCell = 0.5*(xCoordFace1 + xCoordFace2);
        yCoordCell = 0.5*(yCoordFace1 + yCoordFace2);

        // Write the coordinates of the cell centroid and its volume in a tecplot file
        file << std::fixed << std::setprecision(6) << volume[i] << std::endl;
    }

    // Close the tecplot file
    file.close();

    // Inform the user that the file has been created
    std::cout << "The tecplot file " << filename << " has been created." << std::endl;
}