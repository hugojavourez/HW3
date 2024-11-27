#include "math.h"

#include <string>
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

    // Loop over the cells
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

void BoundaryConditions(const int n, const double MachNumber, const double AoA, double fluidProperties[5], const std::vector<int>& cellType, std::vector<int>& cellToFaces, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R) {
    int physicalCellIndex; // Index of the physical cell touching the boundary
    int faceIndex; // Index of the face between the physical and the ghost cell
    double normalVector[2]; // Normal vector of the face

    double dotProductResult = 0; // Dot product of the farfield flow velocity and the normal vector

    double soundSpeed; // Speed of sound

    double flowVelocity[2]; // Flow velocity of the physical cell or the farfield state (in the else part of the loop)
    double velocityMagnitude; // Magnitude of the flow velocity
    double totalEnergy; // Total energy of the physical cell or the farfield state (in the else part of the loop)
    double pPhysical; // Physical cell pressure
    double pGhost, rhoGhost, uGhost, vGhost, EGhost; // Ghost cell properties

    // Loop over the cells
    for (int c = 0; c < cellType.size(); c++) {
        // Check if the cell is a boundary cell (-1 if it is a wall, 1 if it is a farfield)
        if (cellType[c] == -1) {
            // Determine the index of the physical cell touching the boundary
            if (cellType[c+n-1] == -1) { // The current cell is on the second layer of ghost cells
                physicalCellIndex = c + (n - 1) * 2;
            } else { // The current cell is on the first layer of ghost cells
                physicalCellIndex = c + (n - 1);
            }

            // Determine the index of the face between the physical and the ghost cell
            faceIndex = cellToFaces[4*physicalCellIndex];

            // Calculate the flow velocity in the physical cell
            flowVelocity[0] = W[3 * physicalCellIndex + 1] / W[3 * physicalCellIndex];
            flowVelocity[1] = W[3 * physicalCellIndex + 2] / W[3 * physicalCellIndex];
            velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));
            // Calculate the flow velocity in the ghost cell
            uGhost = flowVelocity[0] - 2*velocityMagnitude * xNormal[faceIndex]; 
            vGhost = flowVelocity[1] - 2*velocityMagnitude * yNormal[faceIndex];
            // Calculate the pressure and the total energy in the physical cell
            pGhost = (fluidProperties[2] - 1) * ((W[3 * physicalCellIndex + 3] / W[3 * physicalCellIndex]) - 0.5 * W[3 * physicalCellIndex] * pow(velocityMagnitude,2));
            EGhost = pGhost / ((fluidProperties[2] - 1)) + 0.5 * W[3 * physicalCellIndex] * (pow(uGhost,2) + pow(vGhost,2));

            // Set the values of the ghost cell
            W[3 * c] = W[3 * physicalCellIndex];
            W[3 * c + 1] = W[3 * c] * uGhost;
            W[3 * c + 2] = W[3 * c] * vGhost;
            W[3 * c + 3] = W[3 * c] * EGhost;

            // Setting the residual to zero
            R[c] = 0;

        } else if (cellType[c] == 1) {
            // Determine the index of the physical cell touching the boundary
            if (cellType[c-(n-1)] == -1) { // The current cell is on the second layer of ghost cells
                physicalCellIndex = c - (n - 1) * 2;
            } else { // The current cell is on the first layer of ghost cells
                physicalCellIndex = c - (n - 1);
            }

            // Determine the index of the face between the physical and the ghost cell
            faceIndex = cellToFaces[4*physicalCellIndex];

            // Calculate the flow velocity vector and the normal vector of the face
            flowVelocity[0] = MachNumber * cos(AoA);
            flowVelocity[1] = MachNumber * sin(AoA);
            velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));
            normalVector[0] = xNormal[c];
            normalVector[1] = yNormal[c];

            // Determine if the flow is entering or leaving the domain
            dotProduct(flowVelocity, normalVector, dotProductResult);

            // Determine farfield boundary type
            if (MachNumber < 1.0) {
                if (dotProductResult < 0.0) {
                    // Subsonic inflow
                    soundSpeed = sqrt(fluidProperties[2] * fluidProperties[3] * fluidProperties[4]);
                    pPhysical = (fluidProperties[2] - 1) * ((W[3 * physicalCellIndex + 3] / W[3 * physicalCellIndex]) - 0.5 * W[3 * physicalCellIndex] * pow(sqrt(pow(W[3 * physicalCellIndex + 1],2) + pow(W[3 * physicalCellIndex + 2],2)) / W[3 * physicalCellIndex],2));
                    pGhost = 0.5 * (fluidProperties[1] + pPhysical - W[3 * physicalCellIndex] *soundSpeed * ((flowVelocity[0] - (W[3 * physicalCellIndex + 1] / W[3 * physicalCellIndex])) * xNormal[faceIndex] + (flowVelocity[1] - (W[3 * physicalCellIndex + 2] / W[3 * physicalCellIndex])) * yNormal[faceIndex]));
                    rhoGhost = W[3 * physicalCellIndex] + (pGhost - pPhysical) / pow(soundSpeed,2);
                    uGhost = (W[3 * physicalCellIndex + 1] / W[3 * physicalCellIndex]) - xNormal[faceIndex] * (pPhysical - pGhost) / (W[3 * physicalCellIndex] * soundSpeed);
                    vGhost = (W[3 * physicalCellIndex + 2] / W[3 * physicalCellIndex]) - yNormal[faceIndex] * (pPhysical - pGhost) / (W[3 * physicalCellIndex] * soundSpeed);
                    EGhost = pGhost / ((fluidProperties[2] - 1)) + 0.5 * rhoGhost * (pow(uGhost,2) + pow(vGhost,2));
                } else {
                    // Subsonic outflow
                    soundSpeed = sqrt(fluidProperties[2] * fluidProperties[3] * fluidProperties[4]);
                    pPhysical = (fluidProperties[2] - 1) * ((W[3 * physicalCellIndex + 3] / W[3 * physicalCellIndex]) - 0.5 * W[3 * physicalCellIndex] * pow(sqrt(pow(W[3 * physicalCellIndex + 1],2) + pow(W[3 * physicalCellIndex + 2],2)) / W[3 * physicalCellIndex],2));
                    pGhost = fluidProperties[1];
                    rhoGhost = W[3 * physicalCellIndex] + (pGhost - pPhysical) / pow(soundSpeed,2);
                    uGhost = (W[3 * physicalCellIndex + 1] / W[3 * physicalCellIndex]) + xNormal[faceIndex] * (pPhysical - pGhost) / (W[3 * physicalCellIndex] * soundSpeed);
                    vGhost = (W[3 * physicalCellIndex + 2] / W[3 * physicalCellIndex]) + yNormal[faceIndex] * (pPhysical - pGhost) / (W[3 * physicalCellIndex] * soundSpeed);
                    EGhost = pGhost / ((fluidProperties[2] - 1)) + 0.5 * rhoGhost * (pow(uGhost,2) + pow(vGhost,2));
                }
                W[3 * c] = rhoGhost;
                W[3 * c + 1] = rhoGhost * uGhost;
                W[3 * c + 2] = rhoGhost * vGhost;
                W[3 * c + 3] = rhoGhost * EGhost;
            } else {
                if (dotProductResult < 0.0) {
                    // Supersonic inflow (the state is set to the farfield state)
                    totalEnergy = (fluidProperties[1] / ((fluidProperties[2] - 1)) + 0.5 * fluidProperties[0] * pow(velocityMagnitude,2));
                    W[3 * c] = fluidProperties[0];
                    W[3 * c + 1] = fluidProperties[0] * flowVelocity[0];
                    W[3 * c + 2] = fluidProperties[0] * flowVelocity[1];
                    W[3 * c + 3] = fluidProperties[0] * totalEnergy;
                } else {
                    // Supersonic outflow (the state is set to the physical state)
                    W[3 * c] = W[3 * physicalCellIndex];
                    W[3 * c + 1] = W[3 * physicalCellIndex + 1];
                    W[3 * c + 2] = W[3 * physicalCellIndex + 2];
                    W[3 * c + 3] = W[3 * physicalCellIndex + 3];
                }
            }

            // Setting the residual to zero
            R[c] = 0;
        }
    }
}

void Calculateflux(int n, double MachNumber, double AoA, double fluidProperties[5], int faceNumber, int cellNumber, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, length, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& Fc, std::vector<double>& R,std::vector<double>& F, std::vector<double>& G){
    // Il faut tout d'abord itérer sur les arêtes.

    for (int i = 0; i < faceNumber; i++ ){
        // Get les cellules concernés
        int id_cell_left = faceToCellsLeft[i];
        int id_cell_right = faceToCellsRight[i];

        // on calcule les proprièté au cellules de gauche
        double rho_left = W[id_cell_left];
        double u_left = W[id_cell_left + 1]/W[id_cell_left];
        double v_left = W[id_cell_left + 2]/W[id_cell_left];
        double E_left = W[id_cell_left + 3]/W[id_cell_left];
        double p_left = (fluidProperties[2]-1)*(W[id_cell_left + 3] - 0.5*(rho_left*(u_left*u_left + v_left*v_left)));

        // On calcule les proprièté au cellules de droite
        double rho_right = W[id_cell_right];
        double u_right = W[id_cell_right + 1]/W[id_cell_right];
        double v_right = W[id_cell_right + 2]/W[id_cell_right];
        double E_right = W[id_cell_right + 3]/W[id_cell_right];
        double p_right = (fluidProperties[2]-1)*(W[id_cell_right + 3] - 0.5*(rho_right*(u_right*u_right + v_right*v_right)));

        // Ensuite, Calculer les propriètés au arêtes
        double rho = 0.5*(rho_left + rho_right);
        double u = 0.5*(u_left + u_right);
        double v = 0.5*(v_left + v_right);
        double E = 0.5*(E_right + E_left);
        double p = 0.5*(p_left + p_right);

        // Maintenant on construit F et G, les flux lié respectivement à x et à y
        F = {
            rho*u, 
            rho *u*u + p,
            rho*u*v,
            u*(rho*(E+p)) };
        G = {
            rho*u, 
            rho *u*v,
            rho*v*v + p,
            v*(rho*(E+p))};
        
        // Finally, we construct the fluxes 
        flux_normal = xNormal[i] *F + yNormal[i] * G

        //lenght = face lenght
        flux = flux_normal * lenght[i]  










    }   

}




    
