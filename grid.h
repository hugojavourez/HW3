#pragma once

#include <iostream>
#include <vector>
#include <sstream>

void readCoordinates(const std::string& filename, const int n, std::vector<double>& xCoords, std::vector<double>& yCoords);
void connectivity(const int n, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<int>& cellToFaces, std::vector<int>& faceToNodes, int& faceNumber, int& cellNumber);
void cellVolume(const int n, const int cellNumber, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& volume);
void faceLength(const int n, const int faceNumber, const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<int>& faceToNodes, std::vector<double>& length);
void faceNormal(const int n, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& xNormal, std::vector<double>& yNormal);
void faceAndCellTypes(const int n, std::vector<int>& faceType, std::vector<int>& cellType);
void showGrid();
