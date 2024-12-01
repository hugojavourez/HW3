#pragma once

#include <iostream>
#include <vector>
#include <sstream>

void volumeTest(const int n, const int cellNumber, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<int>& cellToFaces, std::vector<int>& faceToNodes, std::vector<double>& volume, std::string filename);