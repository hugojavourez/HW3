#include "math.h"
#include "manager.h"

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
 * @param n The grid size (nxn).
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
void Initialization(int n, double MachNumber, double AoA, double fluidProperties[5], int cellNumber, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R) {
    //std::cout << "Initialization function called" << std::endl;

    // Rezise vectors to hold the values
    W.resize(4 * cellNumber);
    R.resize(4 * cellNumber);

    // Calculate the flow velocity in each direction
    double flowDirection[2] = {cos(AoA), sin(AoA)};
    double flowVelocity[2] = {flowDirection[0] * MachNumber, flowDirection[1] * MachNumber};
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
            flowVelocity[0] = W[4 * physicalCellIndex + 1] / W[4 * physicalCellIndex];
            flowVelocity[1] = W[4 * physicalCellIndex + 2] / W[4 * physicalCellIndex];
            velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));
            // Calculate the flow velocity in the ghost cell
            uGhost = flowVelocity[0] - 2*velocityMagnitude * xNormal[faceIndex]; 
            vGhost = flowVelocity[1] - 2*velocityMagnitude * yNormal[faceIndex];
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
            if (cellType[c-(n-1)] == -1) { // The current cell is on the second layer of ghost cells
                physicalCellIndex = c - (n - 1) * 2;
            } else { // The current cell is on the first layer of ghost cells
                physicalCellIndex = c - (n - 1);
            }

            // Determine the index of the face between the physical and the ghost cell
            faceIndex = cellToFaces[4*physicalCellIndex + 3];

            // Calculate the flow velocity vector and the normal vector of the face
            flowVelocity[0] = MachNumber * cos(AoA);
            flowVelocity[1] = MachNumber * sin(AoA);
            velocityMagnitude = sqrt(pow(flowVelocity[0],2) + pow(flowVelocity[1],2));
            normalVector[0] = xNormal[faceIndex];
            normalVector[1] = yNormal[faceIndex];

            // Determine if the flow is entering or leaving the domain
            dotProduct(flowVelocity, normalVector, dotProductResult);

            // Determine farfield boundary type
            if (MachNumber < 1.0) {
                if (dotProductResult < 0.0) {
                    // Subsonic inflow
                    soundSpeed = sqrt(fluidProperties[2] * fluidProperties[3] * fluidProperties[4]);
                    pPhysical = (fluidProperties[2] - 1) * ((W[4 * physicalCellIndex + 3] / W[4 * physicalCellIndex]) - 0.5 * W[4 * physicalCellIndex] * pow(sqrt(pow(W[4 * physicalCellIndex + 1],2) + pow(W[4 * physicalCellIndex + 2],2)) / W[4 * physicalCellIndex],2));
                    pGhost = 0.5 * (fluidProperties[1] + pPhysical - W[4 * physicalCellIndex] *soundSpeed * ((flowVelocity[0] - (W[4 * physicalCellIndex + 1] / W[4 * physicalCellIndex])) * xNormal[faceIndex] + (flowVelocity[1] - (W[4 * physicalCellIndex + 2] / W[4 * physicalCellIndex])) * yNormal[faceIndex]));
                    rhoGhost = W[4 * physicalCellIndex] + (pGhost - pPhysical) / pow(soundSpeed,2);
                    uGhost = (W[4 * physicalCellIndex + 1] / W[4 * physicalCellIndex]) - xNormal[faceIndex] * (pPhysical - pGhost) / (W[4 * physicalCellIndex] * soundSpeed);
                    vGhost = (W[4 * physicalCellIndex + 2] / W[4 * physicalCellIndex]) - yNormal[faceIndex] * (pPhysical - pGhost) / (W[4 * physicalCellIndex] * soundSpeed);
                    EGhost = pGhost / ((fluidProperties[2] - 1)) + 0.5 * rhoGhost * (pow(uGhost,2) + pow(vGhost,2));
                } else {
                    // Subsonic outflow
                    soundSpeed = sqrt(fluidProperties[2] * fluidProperties[3] * fluidProperties[4]);
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
 * @param n The grid size (nxn).
 * @param MachNumber The Mach number of the flow.
 * @param AoA The angle of attack in degrees.
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
    // Il faut tout d'abord itérer sur les arêtes.
    // Déclaration de variables
    double rho_left, u_left, v_left, E_left, p_left;
    double rho_right, u_right, v_right, E_right, p_right;
    double rho, u, v, E, p;
    std::vector<double> F1_Roe(4),F5_Roe(4),F234_Roe(4),A_Roe(4); 
    double F1_Roe_term,F5_Roe_term ;
    double rho_Roe, u_Roe, v_Roe, H_Roe, c_Roe, V_Roe, q2_Roe, H_left, H_right,Delta_V,V_c;

    for (int i = 0; i < faceNumber; i++ ){
        // Check if the face is on a side of the domain
        if (faceToCellsRight[i] == -1 || faceToCellsLeft[i] == -1){
            continue;
        }

        // Get les cellules concernés
        int id_cell_left = faceToCellsLeft[i];
        int id_cell_right = faceToCellsRight[i];

        // on calcule les proprièté au cellules de gauche
        rho_left = W[4*id_cell_left];
        u_left = W[4*id_cell_left + 1]/W[4*id_cell_left];
        v_left = W[4*id_cell_left + 2]/W[4*id_cell_left];
        E_left = W[4*id_cell_left + 3]/W[4*id_cell_left];
        p_left = (fluidProperties[2]-1)*(W[4*id_cell_left + 3] - 0.5*(rho_left*(u_left*u_left + v_left*v_left)));

        // Print the values

        //std::cin.get();

        //On calcule les proprièté au cellules de droite
        rho_right = W[4*id_cell_right];
        u_right = W[4*id_cell_right + 1]/W[4*id_cell_right];
        v_right = W[4*id_cell_right + 2]/W[4*id_cell_right];
        E_right = W[4*id_cell_right + 3]/W[4*id_cell_right];
        p_right = (fluidProperties[2]-1)*(W[4*id_cell_right + 3] - 0.5*(rho_right*(u_right*u_right + v_right*v_right)));

        // Ensuite, Calculer les propriètés à l'arète
        rho = 0.5*(rho_left + rho_right);
        u = 0.5*(u_left + u_right);
        v = 0.5*(v_left + v_right);
        E = 0.5*(E_right + E_left);
        p = 0.5*(p_left + p_right);
        

        // Print the values
        //std::cin.get();
        //std::cout << "rho: " << rho << std::endl;
        //std::cout << "u: " << u << std::endl;
        //std::cout << "v: " << v << std::endl;
        //std::cout << "E: " << E << std::endl;
        //std::cout << "p: " << p << std::endl;

        V_c = u*xNormal[i] + v*yNormal[i]; //produit scalaire

        // ADDING ROE CALCULATION
        // Now we calculate the the properties of Roe
        rho_Roe = std::sqrt(rho_left*rho_right);
        u_Roe = (u_left*std::sqrt(rho_left) + u_right*std::sqrt(rho_right)) / (std::sqrt(rho_left) + std::sqrt(rho_right));
        v_Roe = (v_left*std::sqrt(rho_left) + v_right*std::sqrt(rho_right)) / (std::sqrt(rho_left) + std::sqrt(rho_right));
        H_left = (rho_left*E_left + p_left)/rho_left;
        H_right = (rho_right*E_right + p_right)/rho_right;
        H_Roe = (H_left*std::sqrt(rho_left) + H_right*std::sqrt(rho_right)) / (std::sqrt(rho_left) + std::sqrt(rho_right));
        q2_Roe = u_Roe*u_Roe + v_Roe*v_Roe;
        V_Roe = u_Roe * xNormal[i] + v_Roe * yNormal[i];
        c_Roe = std::sqrt((fluidProperties[2]-1)*(H_Roe-q2_Roe/2));
        Delta_V = (u_right - u_left)*xNormal[i] + (v_right - v_left)*yNormal[i];

        // Now we Build the F1
        F1_Roe_term = std::abs(V_Roe - c_Roe) * (((p_right-p_left) - rho_Roe*c_Roe*Delta_V)/(2*c_Roe*c_Roe));
        F1_Roe[0] =  F1_Roe_term;
        F1_Roe[1] =  F1_Roe_term*(u_Roe - c_Roe*xNormal[i]);
        F1_Roe[2] =  F1_Roe_term*(v_Roe - c_Roe*yNormal[i]);
        F1_Roe[3] =  F1_Roe_term*(H_Roe - c_Roe*V_Roe);

        //Now we Build F234
        F234_Roe[0] = V_Roe*(((rho_right-rho_left)- (p_right-p_left)/(c_Roe*c_Roe))*1 + rho_Roe*0);
        F234_Roe[1] = V_Roe*(((rho_right-rho_left)- (p_right-p_left)/(c_Roe*c_Roe))*u_Roe + rho_Roe*((u_right - u_left) - (Delta_V*xNormal[i])));
        F234_Roe[2] = V_Roe*(((rho_right-rho_left)- (p_right-p_left)/(c_Roe*c_Roe))*v_Roe + rho_Roe*((v_right - v_left) - (Delta_V*yNormal[i])));
        F234_Roe[3] = V_Roe*(((rho_right-rho_left)- (p_right-p_left)/(c_Roe*c_Roe))*(q2_Roe/2) + rho_Roe*(u_Roe*(u_right-u_left) + v_Roe*(v_right-v_left) - V_Roe*Delta_V));

        //Now we Build F5
        F5_Roe_term = std::abs(V_Roe + c_Roe) * (((p_right+p_left) - rho_Roe*c_Roe*Delta_V)/(2*c_Roe*c_Roe));
        F5_Roe[0] =  F5_Roe_term;
        F5_Roe[1] =  F5_Roe_term*(u_Roe +  c_Roe*xNormal[i]);
        F5_Roe[2] =  F5_Roe_term*(v_Roe + c_Roe*yNormal[i]);
        F5_Roe[3] =  F5_Roe_term*(H_Roe + c_Roe*V_Roe);


        // Now we build the A_Roe
        for (int i = 0; i<A_Roe.size();i++){
            A_Roe[i] = F1_Roe[i] + F234_Roe[i] + F5_Roe[i];
        }

        // Maintenant on construit F et G, les flux lié respectivement à x et à y

        //std::vector<double> F = {
        //xNormal[i], 
        //yNormal[i],
        //0,
        //0};


        std::vector<double> F = {
            rho*V_c, 
            rho *u*V_c + p*xNormal[i],
            rho*v*V_c + p*yNormal[i],
            u*(rho*(E+p)) };
        // We now substrat the row term, 
        for (int i = 0; i<A_Roe.size();i++){
            F[i] -= 0.5*A_Roe[i];
        }


        // Finally, we construct the fluxes Flux_normal = xNormal[i] *F + yNormal[i] * G;


        for (double &Element : F){
            Element *= length[i]; // on multiplie flux par face lenght 
        }

        // now we must calculate the Residu (LC) and add / remove from the face
        int index_res = 0;
        for (double &Element : F){
            double Left_residual = Element/cellvolume[id_cell_left];
            double Right_residual = Element/cellvolume[id_cell_right];
            R[4*id_cell_left+index_res] -= Left_residual;
            R[4*id_cell_right+index_res] += Right_residual;
            index_res++;
        };
        
        

    }   
} 

void Euler(double dt ,double t, int n, double MachNumber, double AoA, double fluidProperties[5], std::vector<int>& celltype, std::vector<int>& cellToFaces, int faceNumber, int cellNumber, std::vector<double>& cellvolume, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& length, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R){
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

    //std::vector<double> rho(cellNumber);
    //std::vector<double> u(cellNumber);
    //std::vector<double> v(cellNumber);
    //std::vector<double> E(cellNumber);
    //std::vector<double> p(cellNumber);
    
    while (live_time < t){
        // NOTE THAT R HAS ALREADY BEEN DIVIDED BY THE CELL AREA AND MULTIPLIED BY FACE LENGHT
        W_0 = W;
        
        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_0 = R;
        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i]-(dt)*R_0[i]; //W(1)
        }
        live_time += dt;


        // Now that we have W(n+1), lets calculate our properties
        //for (int i = 0 ; i<cellNumber; i++){
            //rho[i] = W[4*i];
            //u[i] = W[4*i + 1] / W[4*i];
            //v[i] = W[4*i + 2] / W[4*i];
            //E[i] = W[4*i + 3] / W[4*i];
            //p[i] = (fluidProperties[2]-1)*(E[i] - 0.5*(rho[i]*(u[i]*u[i] + v[i]*v[i])));

        // Print the residuals R
        //for (int i = 0; i < cellNumber; i++) {
            //std::cout << "R[" << i << "] = " << R[i] << std::endl;
        //}

        convergenceManager(k, cellNumber, R, globalResidual);
        k += 1;
    }
    std::cout << "End of the simulation" << std::endl;
}




void RK4(double dt ,double t, int n, double MachNumber, double AoA, double fluidProperties[5], std::vector<int>& celltype, std::vector<int>& cellToFaces, int faceNumber, int cellNumber, std::vector<double>& cellvolume, std::vector<int>& faceToCellsLeft, std::vector<int>& faceToCellsRight, std::vector<double>& length, std::vector<double>& xNormal, std::vector<double>& yNormal, std::vector<double>& W, std::vector<double>& R){
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

    //std::vector<double> rho(cellNumber);
    //std::vector<double> u(cellNumber);
    //std::vector<double> v(cellNumber);
    //std::vector<double> E(cellNumber);
    //std::vector<double> p(cellNumber);
    
    while (live_time < t){
        // NOTE THAT R HAS ALREADY BEEN DIVIDED BY THE CELL AREA AND MULTIPLIED BY FACE LENGHT
        W_0 = W;
        
        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_0 = R;
        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i]-(dt/2)*R_0[i]; //W(1)
        }
        
        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_1 = R;

        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i]-(dt/2)*R_1[i]; // W(2)
        }
        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_2 = R;
        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i]-(dt)*R_2[i]; // W(3)
        }

        CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, R);
        BoundaryConditions(n, MachNumber, AoA, fluidProperties,celltype, cellToFaces, xNormal, yNormal, W, R);
        R_3 = R;
        for (int i = 0 ; i<cellNumber*4; i++){
            W[i] = W_0[i] - (dt/6) * (R_0[i] + R_1[i] + R_2[i] + R_3[i]); //W(4) = W(n+1)
        }
        live_time += dt;


        // Now that we have W(n+1), lets calculate our properties
        //for (int i = 0 ; i<cellNumber; i++){
            //rho[i] = W[4*i];
            //u[i] = W[4*i + 1] / W[4*i];
            //v[i] = W[4*i + 2] / W[4*i];
            //E[i] = W[4*i + 3] / W[4*i];
            //p[i] = (fluidProperties[2]-1)*(E[i] - 0.5*(rho[i]*(u[i]*u[i] + v[i]*v[i])));

        // Print the residuals R
        //for (int i = 0; i < cellNumber; i++) {
            //std::cout << "R[" << i << "] = " << R[i] << std::endl;
        //}

        convergenceManager(k, cellNumber, R, globalResidual);
        k += 1;
    }
    std::cout << "End of the simulation" << std::endl;
}

/**
 * Writes the Tecplot file for the variable.
 * 
 * @param filename The name of the file to write.
 * @param NI
 * @param NJ
 * @param X The x-coordinates of the grid.
 * @param Y The y-coordinates of the grid.
 * @param variable The variable to write.
 */
void WriteTecplotFile(const std::string & filename,
                      int NI, int NJ,
                      const std::vector<double>& X,
                      const std::vector<double>& Y,
                      const std::vector<double>& variable){
    
    int NCI = NI - 1;
    int NCJ = NJ - 1;

    std::ofstream outfile(filename);

    // Write header
    outfile << "TITLE = \"CFD Simulation Results\"\n";
    outfile << "VARIABLES = \"X\", \"Y\", \"Variable\" \n";
    outfile << "ZONE T=\"Flow Field\", I=" << NI << ", J=" << NJ << ", DATAPACKING=BLOCK\n";
    outfile << "VARLOCATION=([3-7]=CELLCENTERED)\n";

    // Crite the coordinates
    for (int idx = 0; idx < NI * NJ; ++idx) {
        outfile << X[idx] << " ";
    }
    outfile << "\n";

    for (int idx = 0; idx < NI * NJ; ++idx) {
        outfile << Y[idx] << " ";
    }
    outfile << "\n";

    // Write the cell variable
    auto writeCellVariable = [&](const std::vector<double>& var) {
        for (int idx = 0; idx < NCI * NCJ; ++idx) {
            outfile << var[idx] << " ";
        }
        outfile << "\n";
    };
    writeCellVariable(variable);

    outfile.close();
}
