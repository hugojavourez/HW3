#include "math.h"

#include <iostream>

void crossProduct(const double A[2], const double B[2], double result){
    // Cross product calculation
    result = A[0] * B[1] - A[1] * B[0]; // Ax * By - Ay * Bx
}

void norm(const double A[2], double result) {
    // Norm calculation
    result = std::sqrt(std::pow(A[0],2) + std::pow(A[1],2)); // sqrt(Ax² + Ay²)
}