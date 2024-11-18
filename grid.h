#pragma once

#include <iostream>
#include <vector>
#include <sstream>

void readCoordinates(const std::string& filename, std::vector<double>& xCoords, std::vector<double>& yCoords, int n1, int n2);
void cellsVolume(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& volume);
void faceLength(int n1, int n2, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& length);
void showGrid();
