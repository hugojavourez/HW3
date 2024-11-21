#pragma once

#include <iostream>
#include <vector>
#include <sstream>

void readCoordinates(const std::string& filename, std::vector<double>& xCoords, std::vector<double>& yCoords, int n1, int n2);
void cellVolume(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& volume);
void faceLength(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& length);
void faceNormal(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& xNormal, std::vector<double>& yNormal);
void connectivity(int n1, int n2, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, int faceNumber, int cellNumber);
void showGrid();
