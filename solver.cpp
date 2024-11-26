#include "math.h"

#include <iostream>
#include <vector>

void Initialization(int n, double MachNumber, double AoA, double fluidProperties[5], int cellNumber, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R) {
    std::cout << "Initialization function called" << std::endl;

    // Rezise vectors to hold the values
    W.resize(4 * cellNumber);
    R.resize(cellNumber);

    // Calculate the flow velocity in each direction
    double flowDirection[2] = {cos(AoA), sin(AoA)};
    double flowVelocity[2] = {flowDirection[0] * MachNumber, flowDirection[1] * MachNumber};
    double velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));

    // Calculate the total energy
    double totalEnergy = (fluidProperties[1] / ((fluidProperties[2] - 1)) + 0.5 * fluidProperties[0] * pow(velocityMagnitude,2));

    // Calculate all vectors
    for (int c = 0; c < cellNumber; c++) {
        // Calculate the fluid properties vector W
        W[3 * c] = fluidProperties[0];
        W[3 * c + 1] = fluidProperties[0] * flowVelocity[0];
        W[3 * c + 2] = fluidProperties[0] * flowVelocity[1];
        W[3 * c + 3] = fluidProperties[0] * totalEnergy;

        // Calculate the residual vector R (initially set to zero)
        R[c] = 0;
    }
}

void Calculateflux(int n, double MachNumber, double AoA, double fluidProperties[5], int faceNumber, int cellNumber, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& Fc, std::vector<double>& R){
    // Il faut tout d'abord itérer sur les arêtes.

    for (int i = 0; i < faceNumber; i++ ){
        //get les cellules concernés
        int id_cell_left = faceToCellsLeft[i];
        int id_cell_right = faceToCellsRight[i];

        // on calcule les proprièté au cellules de gauche
        double rho_left = W[id_cell_left];
        double u_left = W[id_cell_left + 1]/W[id_cell_left];
        double v_left = W[id_cell_left + 2]/W[id_cell_left];
        double E_left = W[id_cell_left + 3]/W[id_cell_left];
        double p_left = (fluidProperties[2]-1)*(W[id_cell_left + 3] - 0.5*(rho_left*(u_left*u_left + v_left*v_left)));

        // on calcule les proprièté au cellules de droite
        double rho_right = W[id_cell_right];
        double u_right = W[id_cell_right + 1]/W[id_cell_right];
        double v_right = W[id_cell_right + 2]/W[id_cell_right];
        double E_right = W[id_cell_right + 3]/W[id_cell_right];
        double p_right = (fluidProperties[2]-1)*(W[id_cell_right + 3] - 0.5*(rho_right*(u_right*u_right + v_right*v_right)));

        //ensuite, construire la matrice de flux

    }   

}




    
