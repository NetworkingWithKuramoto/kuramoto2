#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "normalizer.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double norm_softsign(double x) {
    return x / (fabs(x) + 1);
}

double norm_fcr(double x, double lower_bound, double upper_bound) {
    double slope, shift;
    slope = 2 / (upper_bound - lower_bound);
    shift = -1 * (upper_bound + lower_bound) / (upper_bound - lower_bound);
    return fmin(1, fmax(-1, slope * x + shift));
}

double norm_linear(double x, double a) {
    if (a == 0.0) {
        // if (x > 0.0) return 1.0;
        // else if (x < 0.0) return -1.0;
        // else return 0.0;
        double eps = 1e-4;
        return x / sqrt(x*x + eps*eps);
    }

    double abs_x = fabs(x);
    return x / ((abs_x > a) ? abs_x : a);
}

double norm_tanh(double x) {
    return tanh(x);
}

double wrap_to_pi(double x) {
    double y = fmod(x + M_PI, 2*M_PI);
    if (y < 0) y += 2*M_PI;
    return fabs(y - M_PI);
}
