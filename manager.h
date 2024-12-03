#pragma once

#include <iostream>
#include <vector>
#include <sstream>

void timeManager(const std::string& caseName, double& startTime, double& endTime);
void convergenceManager(const int iterationNumber, const int cellNumber, const std::vector<double>& R, std::vector<double>& globalResidual);
