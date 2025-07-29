#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITER 10
#define SAFETY 0.8        // Reduced safety factor to avoid aggressive step size growth
#define EPSILON 1.0e-18    // Stricter error tolerance
#define MIN_STEP 1.0e-5  // Minimum step size
#define MAX_STEP 1.0e-3      // Prevents excessively large steps

typedef void (*ODEFunc)(double t, double y[], double yp[]);

/**
 * @brief Allocate a 2D matrix of doubles.
 * @param rows Number of rows.
 * @param cols Number of columns.
 * @return Pointer to allocated matrix.
 */
double **allocate_mat(int rows, int cols) {
    double **matrix = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double *)malloc(cols * sizeof(double));
    }
    return matrix;
}

/**
 * @brief Free a dynamically allocated 2D matrix.
 * @param matrix Pointer to the matrix.
 * @param rows Number of rows.
 */
void free_mat(double **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

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
void modified_midpoint(ODEFunc f, double t, double y[], double h, int neqn, int k, double y_out[]) {
    int i, j, m = 2 * k;
    double **ym = allocate_mat(m + 1, neqn);
    double *yp = (double *)malloc(neqn * sizeof(double));

    f(t, y, yp);
    for (i = 0; i < neqn; i++)
        ym[0][i] = y[i];

    double h2 = h / m;
    for (i = 0; i < neqn; i++)
        ym[1][i] = y[i] + h2 * yp[i];

    for (j = 1; j < m; j++) {
        f(t + j * h2, ym[j], yp);
        for (i = 0; i < neqn; i++)
            ym[j + 1][i] = ym[j - 1][i] + 2.0 * h2 * yp[i];
    }

    f(t + h, ym[m], yp);
    for (i = 0; i < neqn; i++)
        y_out[i] = 0.5 * (ym[m][i] + ym[m - 1][i] + h2 * yp[i]);

    free(yp);
    free_mat(ym, m + 1);
}

/**
 * @brief Uses Richardson Extrapolation to improve accuracy.
 * @param table Table of results from midpoint method.
 * @param k Current extrapolation step.
 * @param y_out Output vector after extrapolation.
 * @param neqn Number of equations.
 */
void extrapolate(double table[MAX_ITER][MAX_ITER], int k, double y_out[], int neqn) {
    for (int i = 1; i <= k; i++) {
        double factor = (double)(4 * i * i) / (4 * i * i - 1);
        for (int j = 0; j < neqn; j++)
            table[k][j] = factor * table[k][j] - (factor - 1) * table[k - 1][j];
    }
    for (int i = 0; i < neqn; i++)
        y_out[i] = table[k][i];
}

/**
 * @brief Gragg-Bulirsch-Stoer Integrator compatible with r8_rkf45.
 * 
 * @param f ODE function of the form f(t, y[], yp[]).
 * @param neqn Number of equations.
 * @param y Solution vector (input/output).
 * @param yp Derivative vector (output).
 * @param t Current time (input/output).
 * @param tout Target time.
 * @param relerr Relative error tolerance.
 * @param abserr Absolute error tolerance.
 * @param flag Integration control flag.
 * @return 0 on success, -1 on failure.
 */
int gbs_integrator(
    void (*f)(double t, double y[], double yp[]),
    int neqn, double y[], double yp[],
    double *t, double tout, double *relerr, double abserr, int *flag) {

    int k;
    double y_temp[MAX_ITER][MAX_ITER];
    double error, h = fmin(fabs(tout - *t), MAX_STEP); // Use a bounded step size
    double h_next = h;
    int step_successful = 0;

    do {
        for (k = 0; k < MAX_ITER; k++) {
            modified_midpoint(f, *t, y, h_next, neqn, k + 1, y_temp[k]);
            if (k > 0)
                extrapolate(y_temp, k, y, neqn);

            if (k > 0) {
                error = 0.0;
                for (int i = 0; i < neqn; i++)
                    error = fmax(error, fabs(y[i] - y_temp[k - 1][i]));

                if (error < (*relerr) * fabs(y[0]) + abserr) {
                    step_successful = 1;
                    break;
                }
            }
        }

        if (!step_successful) {
            h_next *= SAFETY * pow((*relerr * fabs(y[0]) + abserr) / (error + 1e-12), 0.25);
            h_next = fmax(h_next, MIN_STEP);
        } else {
            *t += h_next;
            f(*t, y, yp);
            if (*t >= tout) {
                *flag = 2;
                return 0;
            }
        }
    } while (!step_successful);

    *flag = 1;
    return 0;
}
