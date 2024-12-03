#include "manager.h"

#include <iostream>
#include <vector>
#include <sstream>

/**
 * Manages the time needed to run the simulation.
 * 
 * @param iterationNumber The number of the current iteration.
 * @param caseName The name of the case to manage: "Start", "End" or "Show".
 * @param startTime The time when the simulation started.
 * @param endTime The time when the simulation ended.
 * 
 * The function will print the time needed to run the simulation when the caseName is "End" or "Show".
 * The function will print a message when the caseName is "Start".
 * The function will print an error message when the caseName is invalid.
 * 
 * Warning: startTime and endTime are using the clock() function from the ctime library. They are not real time values.
 */
void timeManager(const int iterationNumber, const std::string& caseName, double& startTime, double& endTime) {
    if (caseName == "Start") {
        startTime = clock();
        std::cout << "Starting the simulation..." << std::endl;
    } else if (caseName == "End") {
        endTime = clock();
        std::cout << "Simulation finished." << std::endl;
        std::cout << "Time needed to run the simulation: " << (endTime - startTime) / CLOCKS_PER_SEC << " seconds." << std::endl;
    } else if (caseName == "Show") {
        endTime = clock();
        std::cout << "Iteration " << iterationNumber << " - Time since the beggining of the simulation: " << (endTime - startTime) / CLOCKS_PER_SEC << " seconds." << std::endl;
        
    } else {
        std::cout << "Invalid case name." << std::endl;
    }
}

/**
 * Calculates the global residual for the iteration.
 * 
 * @param iterationNumber The number of the current iteration.
 * @param cellNumber The total number of cells.
 * @param R A vector containing the residuals of each cell.
 * @param globalResidual The global residual for the iteration.
 */
void convergenceManager(const int iterationNumber, const int cellNumber, const std::vector<double>& R, std::vector<double>& globalResidual) {
    // Resize the global residual vector
    globalResidual.resize(4,0);
    
    // Loop over the cells to calculate the residuals
    for (int c = 0; c < cellNumber; c++) {
        globalResidual[0] += R[c*4]*R[c*4];
        globalResidual[1] += R[c*4+1]*R[c*4+1];
        globalResidual[2] += R[c*4+2]*R[c*4+2];
        globalResidual[3] += R[c*4+3]*R[c*4+3];
    }

    // Calculate the square root of the residuals
    for (int i = 0; i < 4; i++) {
        globalResidual[i] = sqrt((1 / cellNumber) * globalResidual[i]);
    }

    // Print the global residual
    std::cout << "Iteration " << iterationNumber << " - Residuals: " << std::endl;
    std::cout << "Density: " << globalResidual[0] << std::endl;
    std::cout << "Momentum x: " << globalResidual[1] << std::endl;
    std::cout << "Momentum y: " << globalResidual[2] << std::endl;
    std::cout << "Energy: " << globalResidual[3] << std::endl;
}