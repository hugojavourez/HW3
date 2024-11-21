#include "math.h"

#include <iostream>

/**
 * Calculates the dot product of two 2D vectors.
 * 
 * @param A The first vector.
 * @param B The second vector.
 * @param result The result of the dot product (Ax * Bx + Ay * By).
 */
void dotProduct(const double A[2], const double B[2], double result) {
    // Dot product calculation
    result = A[0] * B[0] + A[1] * B[1]; // Ax * Bx + Ay * By
}

/**
 * Calculates the cross product of two 2D vectors.
 * 
 * @param A The first vector.
 * @param B The second vector.
 * @param result The result of the cross product (Ax * By - Ay * Bx).
 */
void crossProduct(const double A[2], const double B[2], double result){
    // Cross product calculation
    result = A[0] * B[1] - A[1] * B[0]; // Ax * By - Ay * Bx
}

/**
 * Calculates the norm of a 2D vector.
 * 
 * @param A The vector.
 * @param result The result of the norm (sqrt(Ax² + Ay²)).
 */
void norm(const double A[2], double result) {
    // Norm calculation
    result = std::sqrt(std::pow(A[0],2) + std::pow(A[1],2)); // sqrt(Ax² + Ay²)
}