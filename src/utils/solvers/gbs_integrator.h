#ifndef GBS_INTEGRATOR_H
#define GBS_INTEGRATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Constants for numerical stability */
#define MAX_ITER 10
#define SAFETY 0.8        // Safety factor for step size adaptation
#define EPSILON 1.0e-8    // Error tolerance
#define MIN_STEP 1.0e-5  // Minimum step size
#define MAX_STEP 1.0e-1      // Maximum step size

/* Function pointer for ODE system */
typedef void (*ODEFunc)(double t, double y[], double yp[]);

/**
 * @brief Allocates a 2D matrix.
 * @param rows Number of rows.
 * @param cols Number of columns.
 * @return Pointer to allocated matrix.
 */
double **allocate_mat(int rows, int cols);

/**
 * @brief Frees a dynamically allocated 2D matrix.
 * @param matrix Pointer to the matrix.
 * @param rows Number of rows.
 */
void free_mat(double **matrix, int rows);

/**
 * @brief Implements the Modified Midpoint Method for step doubling.
 * @param f Function defining the ODE system.
 * @param t Current time.
 * @param y Current state vector.
 * @param h Step size.
 * @param neqn Number of equations.
 * @param k Current iteration step.
 * @param y_out Output vector after step.
 */
void modified_midpoint(ODEFunc f, double t, double y[], double h, int neqn, int k, double y_out[]);

/**
 * @brief Uses Richardson Extrapolation to improve accuracy.
 * @param table Table of results from midpoint method.
 * @param k Current extrapolation step.
 * @param y_out Output vector after extrapolation.
 * @param neqn Number of equations.
 */
void extrapolate(double table[MAX_ITER][MAX_ITER], int k, double y_out[], int neqn);

/**
 * @brief Gragg-Bulirsch-Stoer Integrator, compatible with r8_rkf45.
 * 
 * This function integrates a system of `neqn` first-order ODEs from `t` to `tout`.
 * It adjusts the step size adaptively to maintain the given relative and absolute error tolerances.
 * 
 * @param f Function pointer defining the system of ODEs.
 * @param neqn Number of equations.
 * @param y State vector (input/output).
 * @param yp Derivative vector (output).
 * @param t Current time (input/output).
 * @param tout Target time.
 * @param relerr Relative error tolerance.
 * @param abserr Absolute error tolerance.
 * @param flag Integration control flag (similar to r8_rkf45).
 * @return 0 on success, -1 on failure (e.g., step size too small).
 */
int gbs_integrator(
    void (*f)(double t, double y[], double yp[]),
    int neqn, double y[], double yp[],
    double *t, double tout, double *relerr, double abserr, int *flag);

#endif /* GBS_INTEGRATOR_H */
