#include "math.h"

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
            W[4 * c + 1] = W[3 * c] * uGhost;
            W[4 * c + 2] = W[3 * c] * vGhost;
            W[4 * c + 3] = W[3 * c] * EGhost;

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
            R[c] = 0;
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
 * @param LC A vector to store the residuals.
 */
void CalculateResidual( double fluidProperties[5], int faceNumber, std::vector<int> &faceToCellsLeft, std::vector<int> &faceToCellsRight, std::vector<double> &length,std::vector<double> &cellvolume, std::vector<double> &xNormal, std::vector<double> &yNormal, std::vector<double> &W, std::vector<double>& LC){
    // Il faut tout d'abord itérer sur les arêtes.

    for (int i = 0; i < faceNumber; i++ ){
        // Get les cellules concernés
        int id_cell_left = faceToCellsLeft[i];
        int id_cell_right = faceToCellsRight[i];

        // on calcule les proprièté au cellules de gauche
        double rho_left = W[4*id_cell_left];
        double u_left = W[4*id_cell_left + 1]/W[4*id_cell_left];
        double v_left = W[4*id_cell_left + 2]/W[4*id_cell_left];
        double E_left = W[4*id_cell_left + 3]/W[4*id_cell_left];
        double p_left = (fluidProperties[2]-1)*(W[4*id_cell_left + 3] - 0.5*(rho_left*(u_left*u_left + v_left*v_left)));

        // On calcule les proprièté au cellules de droite
        double rho_right = W[4*id_cell_right];
        double u_right = W[4*id_cell_right + 1]/W[4*id_cell_right];
        double v_right = W[4*id_cell_right + 2]/W[4*id_cell_right];
        double E_right = W[4*id_cell_right + 3]/W[4*id_cell_right];
        double p_right = (fluidProperties[2]-1)*(W[4*id_cell_right + 3] - 0.5*(rho_right*(u_right*u_right + v_right*v_right)));

        // Ensuite, Calculer les propriètés à l'arète
        double rho = 0.5*(rho_left + rho_right);
        double u = 0.5*(u_left + u_right);
        double v = 0.5*(v_left + v_right);
        double E = 0.5*(E_right + E_left);
        double p = 0.5*(p_left + p_right);

        // Maintenant on construit F et G, les flux lié respectivement à x et à y
        std::vector<double> F = {
            rho*u, 
            rho *u*u + p,
            rho*u*v,
            u*(rho*(E+p)) };
        std::vector<double> G = {
            rho*v, 
            rho *u*v,
            rho*v*v + p,
            v*(rho*(E+p))};
        
        // Finally, we construct the fluxes Flux_normal = xNormal[i] *F + yNormal[i] * G;

        for (double &Element : F){
            Element *= xNormal[i]; //xNormal[i] *F 
        }

        for (double &Element : G){
            Element *= yNormal[i]; // yNormal[i] * G
        }
        
        std::vector<double> Flux_Normal(G.size());

        for (int j = 0; j< Flux_Normal.size();j++){
            Flux_Normal[j] = F[j] + G[j]; //Flux_normal = xNormal[i] *F + yNormal[i] * G
        }


        for (double &Element : Flux_Normal){
            Element *= length[i]; // on multiplie flux par face lenght
            
        }

        // now we must calculate the Residu (LC) and add / remove from the face
        int index_res = 0;
        for (double &Element : Flux_Normal){
            double Left_residual = Element/cellvolume[id_cell_left];
            double Right_residual = Element/cellvolume[id_cell_right];
            LC[4*id_cell_left+index_res] -= Element/cellvolume[id_cell_left];
            LC[4*id_cell_right+index_res] += Element/cellvolume[id_cell_right];
            index_res++;
        };

        //Lc is what we want

    }   
} 

void RK4(double dt,double t,  int n, double MachNumber, double AoA, double fluidProperties[5], int faceNumber, int cellNumber,std::vector<double> &cellvolume, std::vector<int> &faceToCellsLeft, std::vector<int> &faceToCellsRight, std::vector<double> &length, std::vector<double> &xNormal, std::vector<double> &yNormal, std::vector<double> &W,std::vector<double>& LC){
    
    double live_time = 0;
    std::vector<double> W_0;
    std::vector<double> LC_0;
    std::vector<double> LC_1;
    std::vector<double> LC_2;
    std::vector<double> LC_3;
    std::vector<double> rho(cellNumber);
    std::vector<double> u(cellNumber);
    std::vector<double> v(cellNumber);
    std::vector<double> E(cellNumber);
    std::vector<double> p(cellNumber);
    
    while (live_time<0){
        // NOTE THAT LC HAS ALREADY BEEN DIVIDED BY THE CELL AREA AND MULTIPLIED BY FACE LENGHT
        W_0 = W;
        
        for (int i = 0 ; i<cellNumber*4; i++){
            CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, LC); // LC(0)
            LC_0 = LC;
            W[i] = W_0[i]-(dt/2)*LC_0[i]; // W(1)
            CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, LC); // LC(1)
            LC_1 = LC;
            W[i] = W_0[i]-(dt/2)*LC[i]; // W(2)
            CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, LC); // LC(2)
            LC_2 = LC;
            W[i] = W_0[i]-(dt)*LC[i]; // W(3)
            CalculateResidual(fluidProperties,faceNumber, faceToCellsLeft, faceToCellsRight, length, cellvolume, xNormal, yNormal, W, LC); // LC(3)
            LC_3 = LC;
            W[i] = W_0[i] - (dt/6) * (LC_0[i] + LC_1[i] + LC_2[i] + LC_3[i]); //W(4) = W(n+1)
       }

        // Now that we have W(n+1), lets calculate our properties
        for (int i = 0 ; i<cellNumber; i++){
            rho[i] = W[4*i];
            u[i] = W[4*i + 1] / W[4*i];
            v[i] = W[4*i + 2] / W[4*i];
            E[i] = W[4*i + 3] / W[4*i];
            p[i] = (fluidProperties[2]-1)*(E[i] - 0.5*(rho[i]*(u[i]*u[i] + v[i]*v[i])));
        }
        live_time += dt;
    }
}




    
