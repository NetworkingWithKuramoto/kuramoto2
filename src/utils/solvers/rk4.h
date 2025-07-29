#ifndef RK4_H
#define RK4_H

// Function declaration for rk4
void rk4(double y[], double dydx[], int n, double x, double h, double yout[], void (*derivs)(double, double *, double *));

#endif // RK4_H
