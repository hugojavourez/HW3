#include "math.h"
#include "manager.h"
#include "solution.h"

#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>

/**
 * Initializes the flow variables and the residuals.
 * 
 * @param MachNumber The Mach number of the flow.
 * @param AoA The angle of attack in degrees.
 * @param fluidProperties An array containing the fluid properties: [rhoInf, pInf, gamma, R, TInf].
 * @param cellNumber The total number of cells.
 * @param faceToCellsLeft A vector containing the index of the cell on the left of each face.
 * @param faceToCellsRight A vector containing the index of the cell on the right of each face.
 * @param xNormal A vector containing the x-component of the normal vector of each face.
 * @param yNormal A vector containing the y-component of the normal vector of each face.
 * @param W A vector to store the flow variables.
 * @param R A vector to store the residuals.
 */
void Initialization(double MachNumber, double AoA, double fluidProperties[5], int cellNumber, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R) {
    // Rezise vectors to hold the values
    W.resize(4 * cellNumber);
    R.resize(4 * cellNumber);

    // Calculate the flow velocity in each direction
    double flowDirection[2] = {cos(AoA), sin(AoA)};
    double soundSpeed = sqrt(fluidProperties[2] * fluidProperties[3] * fluidProperties[4]);
    double flowVelocity[2] = {flowDirection[0] * MachNumber * soundSpeed, flowDirection[1] * MachNumber * soundSpeed};
    double velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));

    // Calculate the total energy
    double totalEnergy = (fluidProperties[1] / ((fluidProperties[2] - 1)) + 0.5 * fluidProperties[0] * pow(velocityMagnitude,2));

    // Loop over the cells
    for (int c = 0; c < cellNumber; c++) {
        // Calculate the fluid properties vector W
        W[4 * c] = fluidProperties[0];
        W[4 * c + 1] = fluidProperties[0] * flowVelocity[0];
        W[4 * c + 2] = fluidProperties[0] * flowVelocity[1];
        W[4 * c + 3] = fluidProperties[0] * totalEnergy;

        // Calculate the residual vector R (initially set to zero)
        R[4*c] = 0;
        R[4*c + 1] = 0;
        R[4*c + 2] = 0;
        R[4*c + 3] = 0;
    }
}

/**
 * Applies the boundary conditions to the flow variables: farfield and wall.
 * 
 * @param n The grid size (nxn).
 * @param MachNumber The Mach number of the flow.
 * @param AoA The angle of attack in degrees.
 * @param fluidProperties An array containing the fluid properties: [rhoInf, pInf, gamma, R, TInf].
 * @param cellType A vector containing the type of each cell.
 * @param cellToFaces A vector containing the index of the face between the physical and the ghost cell.
 * @param xNormal A vector containing the x-component of the normal vector of each face.
 * @param yNormal A vector containing the y-component of the normal vector of each face.
 * @param W A vector containing the flow variables.
 * @param R A vector containing the residuals.
 * 
 * The boundary conditions are applied to the ghost cells only.
 */
void BoundaryConditions(const int n1, const int n2, const double MachNumber, const double AoA, double fluidProperties[5], const std::vector<int>& cellType, std::vector<int>& cellToFaces, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R) {
    int physicalCellIndex; // Index of the physical cell touching the boundary
    int faceIndex; // Index of the face between the physical and the ghost cell
    double normalVector[2]; // Normal vector of the face

    double dotProductResult = 0; // Dot product of the farfield flow velocity and the normal vector

    double soundSpeed; // Speed of sound

    double flowVelocity[2]; // Flow velocity of the physical cell or the farfield state (in the else part of the loop)
    double velocityMagnitude; // Magnitude of the flow velocity
    double V2n; // Dot product of the flow velocity and the normal vector
    double totalEnergy; // Total energy of the physical cell or the farfield state (in the else part of the loop)
    double pPhysical; // Physical cell pressure
    double pGhost, rhoGhost, uGhost, vGhost, EGhost; // Ghost cell properties

    // Loop over the cells
    for (int c = 0; c < cellType.size(); c++) {
        // Check if the cell is a boundary cell (-1 if it is a wall, 1 if it is a farfield)
        if (cellType[c] == -1) {
            // Determine the index of the physical cell touching the boundary
            if (cellType[c+n1-1] == -1) { // The current cell is on the second layer of ghost cells
                physicalCellIndex = c + (n1 - 1) * 2;
            } else { // The current cell is on the first layer of ghost cells
                physicalCellIndex = c + (n1 - 1);
            }

            // Determine the index of the face between the physical and the ghost cell
            faceIndex = cellToFaces[4*physicalCellIndex];

            // Calculate the flow velocity in the physical cell
            flowVelocity[0] = W[4 * physicalCellIndex + 1] / W[4 * physicalCellIndex];
            flowVelocity[1] = W[4 * physicalCellIndex + 2] / W[4 * physicalCellIndex];
            velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));

            // Calculate the flow velocity in the ghost cell
            normalVector[0] = xNormal[faceIndex];
            normalVector[1] = yNormal[faceIndex];
            dotProduct(flowVelocity, normalVector, V2n);
            uGhost = flowVelocity[0] - 2 * V2n * xNormal[faceIndex]; 
            vGhost = flowVelocity[1] - 2 * V2n * yNormal[faceIndex];

            // Calculate the pressure and the total energy in the physical cell
            pGhost = (fluidProperties[2] - 1) * ((W[4 * physicalCellIndex + 3] / W[4 * physicalCellIndex]) - 0.5 * W[4 * physicalCellIndex] * pow(velocityMagnitude,2));
            EGhost = pGhost / ((fluidProperties[2] - 1)) + 0.5 * W[4 * physicalCellIndex] * (pow(uGhost,2) + pow(vGhost,2));

            // Set the values of the ghost cell
            W[4 * c] = W[4 * physicalCellIndex];
            W[4 * c + 1] = W[4 * c] * uGhost;
            W[4 * c + 2] = W[4 * c] * vGhost;
            W[4 * c + 3] = W[4 * c] * EGhost;

            // Setting the residual to zero
            R[4*c] = 0;
            R[4*c + 1] = 0;
            R[4*c + 2] = 0;
            R[4*c + 3] = 0;

        } else if (cellType[c] == 1) {
            // Determine the index of the physical cell touching the boundary
            if (cellType[c-(n1-1)] == -1) { // The current cell is on the second layer of ghost cells
                physicalCellIndex = c - (n1 - 1) * 2;
            } else { // The current cell is on the first layer of ghost cells
                physicalCellIndex = c - (n1 - 1);
            }

            // Determine the index of the face between the physical and the ghost cell
            faceIndex = cellToFaces[4*physicalCellIndex + 3];

            // Calculate the flow velocity vector and the normal vector of the face
            soundSpeed = sqrt(fluidProperties[2] * fluidProperties[3] * fluidProperties[4]);
            flowVelocity[0] = soundSpeed * MachNumber * cos(AoA);
            flowVelocity[1] = soundSpeed * MachNumber * sin(AoA);
            velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));
            normalVector[0] = xNormal[faceIndex];
            normalVector[1] = yNormal[faceIndex];

            // Determine if the flow is entering or leaving the domain
            dotProduct(flowVelocity, normalVector, dotProductResult);

            // Determine farfield boundary type
            if (MachNumber < 1.0) {
                if (dotProductResult < 0.0) {
                    // Subsonic inflow
                    pPhysical = (fluidProperties[2] - 1) * ((W[4 * physicalCellIndex + 3] / W[4 * physicalCellIndex]) - 0.5 * W[4 * physicalCellIndex] * pow(sqrt(pow(W[4 * physicalCellIndex + 1],2) + pow(W[4 * physicalCellIndex + 2],2)) / W[4 * physicalCellIndex],2));
                    pGhost = 0.5 * (fluidProperties[1] + pPhysical - W[4 * physicalCellIndex] *soundSpeed * ((flowVelocity[0] - (W[4 * physicalCellIndex + 1] / W[4 * physicalCellIndex])) * xNormal[faceIndex] + (flowVelocity[1] - (W[4 * physicalCellIndex + 2] / W[4 * physicalCellIndex])) * yNormal[faceIndex]));
                    rhoGhost = W[4 * physicalCellIndex] + (pGhost - pPhysical) / pow(soundSpeed,2);
                    uGhost = (W[4 * physicalCellIndex + 1] / W[4 * physicalCellIndex]) - xNormal[faceIndex] * (pPhysical - pGhost) / (W[4 * physicalCellIndex] * soundSpeed);
                    vGhost = (W[4 * physicalCellIndex + 2] / W[4 * physicalCellIndex]) - yNormal[faceIndex] * (pPhysical - pGhost) / (W[4 * physicalCellIndex] * soundSpeed);
                    EGhost = pGhost / ((fluidProperties[2] - 1)) + 0.5 * rhoGhost * (pow(uGhost,2) + pow(vGhost,2));
                } else {
                    // Subsonic outflow
                    pPhysical = (fluidProperties[2] - 1) * ((W[4 * physicalCellIndex + 3] / W[4 * physicalCellIndex]) - 0.5 * W[4 * physicalCellIndex] * pow(sqrt(pow(W[4 * physicalCellIndex + 1],2) + pow(W[4 * physicalCellIndex + 2],2)) / W[4 * physicalCellIndex],2));
                    pGhost = fluidProperties[1];
                    rhoGhost = W[4 * physicalCellIndex] + (pGhost - pPhysical) / pow(soundSpeed,2);
                    uGhost = (W[4 * physicalCellIndex + 1] / W[4 * physicalCellIndex]) + xNormal[faceIndex] * (pPhysical - pGhost) / (W[4 * physicalCellIndex] * soundSpeed);
                    vGhost = (W[4 * physicalCellIndex + 2] / W[4 * physicalCellIndex]) + yNormal[faceIndex] * (pPhysical - pGhost) / (W[4 * physicalCellIndex] * soundSpeed);
                    EGhost = pGhost / ((fluidProperties[2] - 1)) + 0.5 * rhoGhost * (pow(uGhost,2) + pow(vGhost,2));
                }
                W[4 * c] = rhoGhost;
                W[4 * c + 1] = rhoGhost * uGhost;
                W[4 * c + 2] = rhoGhost * vGhost;
                W[4 * c + 3] = rhoGhost * EGhost;
            } else {
                if (dotProductResult < 0.0) {
                    // Supersonic inflow (the state is set to the farfield state)
                    totalEnergy = (fluidProperties[1] / ((fluidProperties[2] - 1)) + 0.5 * fluidProperties[0] * pow(velocityMagnitude,2));
                    W[4 * c] = fluidProperties[0];
                    W[4 * c + 1] = fluidProperties[0] * flowVelocity[0];
                    W[4 * c + 2] = fluidProperties[0] * flowVelocity[1];
                    W[4 * c + 3] = fluidProperties[0] * totalEnergy;
                } else {
                    // Supersonic outflow (the state is set to the physical state)
                    W[4 * c] = W[4 * physicalCellIndex];
                    W[4 * c + 1] = W[4 * physicalCellIndex + 1];
                    W[4 * c + 2] = W[4 * physicalCellIndex + 2];
                    W[4 * c + 3] = W[4 * physicalCellIndex + 3];
                }
            }

            // Setting the residual to zero
            R[4*c] = 0;
            R[4*c + 1] = 0;
            R[4*c + 2] = 0;
            R[4*c + 3] = 0;
        }
    }
}

/**
 * Calculates the fluxes.
 * 
 * @param fluidProperties An array containing the fluid properties: [rhoInf, pInf, gamma, R, TInf].
 * @param faceNumber The total number of faces.
 * @param cellNumber The total number of cells.
 * @param faceToCellsLeft A vector containing the index of the cell on the left of each face.
 * @param faceToCellsRight A vector containing the index of the cell on the right of each face.
 * @param xNormal A vector containing the x-component of the normal vector of each face.
 * @param yNormal A vector containing the y-component of the normal vector of each face.
 * @param W A vector containing the flow variables.
 * @param R A vector containing the residuals.
 */
void CalculateResidual( double fluidProperties[5], int faceNumber, std::vector<int> &faceToCellsLeft, std::vector<int> &faceToCellsRight, std::vector<double> &length,std::vector<double> &cellvolume, std::vector<double> &xNormal, std::vector<double> &yNormal, std::vector<double> &W,std::vector<double>& R){
    // Variables declaration
    double rho_left, u_left, v_left, E_left, p_left; // Left cell properties
    double rho_right, u_right, v_right, E_right, p_right; // Right cell properties
    double rho, u, v, E, p; // Roe properties
    std::vector<double> F1_Roe(4),F5_Roe(4),F234_Roe(4),A_Roe(4); // Roe fluxes
    double F1_Roe_term,F5_Roe_term; // Roe fluxes terms
    double rho_Roe, u_Roe, v_Roe, H_Roe, c_Roe, V_Roe, q2_Roe, H_left, H_right,Delta_V,V_c_left,V_c_right; // Roe properties
    double LAMBDA_c; // One eigenvalue of the Jacobian matrix
    double Flux[4]; // Flux vector

    // Calculate the speed of sound
    double speedOfSound = sqrt(fluidProperties[2] * fluidProperties[3] * fluidProperties[4]);

    for (int i = 0; i < faceNumber; i++ ){
        // Check if the face is on a side of the domain
        if (faceToCellsRight[i] == -1 || faceToCellsLeft[i] == -1){
            continue;
        }

        // Get the indexs of the cells on the left and right of the face
        int id_cell_left = faceToCellsLeft[i];
        int id_cell_right = faceToCellsRight[i];

        // Calculate the properties of the left cell
        rho_left = W[4*id_cell_left];
        u_left = W[4*id_cell_left + 1]/W[4*id_cell_left];
        v_left = W[4*id_cell_left + 2]/W[4*id_cell_left];
        E_left = W[4*id_cell_left + 3]/W[4*id_cell_left];
        p_left = (fluidProperties[2]-1)*(W[4*id_cell_left + 3] - 0.5*(rho_left*(u_left*u_left + v_left*v_left)));
        H_left = (rho_left*E_left + p_left)/rho_left;
        V_c_left = u_left*xNormal[i] + v_left*yNormal[i];

        // Calculate the properties of the right cell
        rho_right = W[4*id_cell_right];
        u_right = W[4*id_cell_right + 1]/W[4*id_cell_right];
        v_right = W[4*id_cell_right + 2]/W[4*id_cell_right];
        E_right = W[4*id_cell_right + 3]/W[4*id_cell_right];
        p_right = (fluidProperties[2]-1)*(W[4*id_cell_right + 3] - 0.5*(rho_right*(u_right*u_right + v_right*v_right)));
        H_right = (rho_right*E_right + p_right)/rho_right;
        V_c_right = u_right*xNormal[i] + v_right*yNormal[i];

        // Calculate the properties on the face
        rho = 0.5*(rho_left + rho_right);
        u = 0.5*(u_left + u_right);
        v = 0.5*(v_left + v_right);
        E = 0.5*(E_right + E_left);
        p = 0.5*(p_left + p_right);
        Delta_V = (u_right - u_left)*xNormal[i] + (v_right - v_left)*yNormal[i];

        // ROE CALCULATION
        // Calculate the properties of Roe
        rho_Roe = sqrt(rho_left*rho_right);
        u_Roe = (u_left*sqrt(rho_left) + u_right*sqrt(rho_right)) / (sqrt(rho_left) + sqrt(rho_right));
        v_Roe = (v_left*sqrt(rho_left) + v_right*sqrt(rho_right)) / (sqrt(rho_left) + sqrt(rho_right));
        H_Roe = (H_left*sqrt(rho_left) + H_right*sqrt(rho_right)) / (sqrt(rho_left) + sqrt(rho_right));
        q2_Roe = u_Roe*u_Roe + v_Roe*v_Roe;
        V_Roe = u_Roe * xNormal[i] + v_Roe * yNormal[i];
        c_Roe = sqrt((fluidProperties[2]-1)*(H_Roe-q2_Roe/2));

        // Build the F1 Roe term
        if (abs(V_Roe - c_Roe) <= speedOfSound/10) { // Entropy correction
            LAMBDA_c = ((V_Roe - c_Roe)*(V_Roe - c_Roe) + (speedOfSound/10)*(speedOfSound/10)) / (2*(speedOfSound/10));
        } else {
            LAMBDA_c = abs(V_Roe - c_Roe);
        }
        F1_Roe_term = LAMBDA_c * (((p_right-p_left) - rho_Roe*c_Roe*Delta_V)/(2*c_Roe*c_Roe));
        F1_Roe[0] =  F1_Roe_term;
        F1_Roe[1] =  F1_Roe_term*(u_Roe - c_Roe*xNormal[i]);
        F1_Roe[2] =  F1_Roe_term*(v_Roe - c_Roe*yNormal[i]);
        F1_Roe[3] =  F1_Roe_term*(H_Roe - c_Roe*V_Roe);

        // Build the F234 Roe terms
        if (V_Roe <= speedOfSound/10) { // Entropy correction
            LAMBDA_c = (V_Roe*V_Roe + (speedOfSound/10)*(speedOfSound/10)) / (2*(speedOfSound/10));
        } else {
            LAMBDA_c = abs(V_Roe);
        }
        F234_Roe[0] = LAMBDA_c*(((rho_right-rho_left)- (p_right-p_left)/(c_Roe*c_Roe))*1 + rho_Roe*0);
        F234_Roe[1] = LAMBDA_c*(((rho_right-rho_left)- (p_right-p_left)/(c_Roe*c_Roe))*u_Roe + rho_Roe*((u_right - u_left) - (Delta_V*xNormal[i])));
        F234_Roe[2] = LAMBDA_c*(((rho_right-rho_left)- (p_right-p_left)/(c_Roe*c_Roe))*v_Roe + rho_Roe*((v_right - v_left) - (Delta_V*yNormal[i])));
        F234_Roe[3] = LAMBDA_c*(((rho_right-rho_left)- (p_right-p_left)/(c_Roe*c_Roe))*(q2_Roe/2) + rho_Roe*(u_Roe*(u_right-u_left) + v_Roe*(v_right-v_left) - V_Roe*Delta_V));

        // Build of F5 Roe term
        if (abs(V_Roe + c_Roe) <= speedOfSound/10) { // Entropy correction
            LAMBDA_c = ((V_Roe + c_Roe)*(V_Roe + c_Roe) + (speedOfSound/10)*(speedOfSound/10)) / (2*(speedOfSound/10));
        } else {
            LAMBDA_c = abs(V_Roe + c_Roe);
        }
        F5_Roe_term = LAMBDA_c * (((p_right-p_left) + rho_Roe*c_Roe*Delta_V)/(2*c_Roe*c_Roe));
        F5_Roe[0] =  F5_Roe_term;
        F5_Roe[1] =  F5_Roe_term*(u_Roe +  c_Roe*xNormal[i]);
        F5_Roe[2] =  F5_Roe_term*(v_Roe + c_Roe*yNormal[i]);
        F5_Roe[3] =  F5_Roe_term*(H_Roe + c_Roe*V_Roe);

        // Build the A_Roe
        for (int j = 0; j < 4; j++){
            A_Roe[j] = abs(F1_Roe[j]) + abs(F234_Roe[j]) + abs(F5_Roe[j]);
        }

        //std::vector<double> F = {
        //xNormal[i], 
        //yNormal[i],
        //0,
        //0};

        // Calculate F and G, the fluxes related to the left and right cells
        std::vector<double> F = {
            0.5*(rho_left*V_c_left), 
            0.5*(rho_left*u_left*V_c_left + p_left*xNormal[i]),
            0.5*(rho_left*v_left*V_c_left + p_left*yNormal[i]),
            0.5*(rho_left*E_left*V_c_left + p_left*V_c_left)};
        std::vector<double> G = {
            0.5*(rho_right*V_c_right), 
            0.5*(rho_right*u_right*V_c_right + p_right*xNormal[i]),
            0.5*(rho_right*v_right*V_c_right + p_right*yNormal[i]),
            0.5*(rho_right*E_right*V_c_right + p_right*V_c_right)};
        for (int j = 0; j < 4; j++){
            Flux[j] = 0.5*(F[j] + G[j] - A_Roe[j]);
        }

        // Calculate the residual and add / remove from the adjacent cells
        for (int j = 0; j < 4; j++){
            R[4*id_cell_left+j] -= Flux[j]*length[i]/cellvolume[id_cell_left];
            R[4*id_cell_right+j] += Flux[j]*length[i]/cellvolume[id_cell_left];
        };
    }  
} 

void Euler(double dt ,double t, int n1, int n2, double MachNumber, double AoA, double fluidProperties[5], std::vector<int>& celltype, std::vector<int>& cellToFaces, int faceNumber, int cellNumber, std::vector<double>& cellvolume, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& length, std::vector<double>& xCoords, std::vector<double>& yCoords, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R){
    double live_time = 0;
    std::vector<double> W_0;
    std::vector<double> R_0;


    int k = 0; // Iteration number
    double globalResidual = 1; // Residual used for the global convergence

    // Resizing the vectors
    // W_0.resize(4*cellNumber,0.0);
    // R_0.resize(cellNumber,0.0);
    // R_1.resize(cellNumber,0.0);
    // R_2.resize(cellNumber,0.0);
    // R_3.resize(cellNumber,0.0);

    while (live_time < t){
        // NOTE THAT R HAS ALREADY BEEN DIVIDED BY THE CELL AREA AND MULTIPLIED BY FACE LENGHT
        W_0 = W;
        
        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n1, n2, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_0 = R;
        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i]+(dt)*R_0[i]; //W(1)
        }
        live_time += dt;

        std::vector<double> rho, u, v, VMag, E, p;
        writeProperties(n1, n2, xCoords, yCoords, W, fluidProperties, cellNumber, rho, u, v, VMag, E, p);

        convergenceManager(k, cellNumber, R, globalResidual);
        k += 1;
    }
    std::cout << "End of the simulation" << std::endl;
}

void RK4(double dt ,double t, int n1, int n2, double MachNumber, double AoA, double fluidProperties[5], std::vector<int>& celltype, std::vector<int>& cellToFaces, int faceNumber, int cellNumber, std::vector<double>& cellvolume, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& length, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R){
    double live_time = 0;
    std::vector<double> W_0;
    std::vector<double> R_0;
    std::vector<double> R_1;
    std::vector<double> R_2;
    std::vector<double> R_3;

    int k = 0; // Iteration number
    double globalResidual = 1; // Residual used for the global convergence

    // Resizing the vectors
    // W_0.resize(4*cellNumber,0.0);
    // R_0.resize(cellNumber,0.0);
    // R_1.resize(cellNumber,0.0);
    // R_2.resize(cellNumber,0.0);
    // R_3.resize(cellNumber,0.0);
    
    while (live_time < t){
        // NOTE THAT R HAS ALREADY BEEN DIVIDED BY THE CELL AREA AND MULTIPLIED BY FACE LENGHT
        W_0 = W;
        
        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n1, n2, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_0 = R;
        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i]-(dt/2)*R_0[i]; //W(1)
        }
        
        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n1, n2, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_1 = R;

        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i]-(dt/2)*R_1[i]; // W(2)
        }
        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n1, n2, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_2 = R;
        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i]-(dt)*R_2[i]; // W(3)
        }

        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n1, n2, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_3 = R;
        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i] - (dt/6) * (R_0[i] + R_1[i] + R_2[i] + R_3[i]); //W(4) = W(n+1)
        }
        live_time += dt;

        convergenceManager(k, cellNumber, R, globalResidual);
        k += 1;
    }
    std::cout << "End of the simulation" << std::endl;
}
