#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "normalizer.h"

double norm_tanh(double x) {
    return x / (fabs(x) + 1);
}

double norm_fcr(double x, double lower_bound, double upper_bound) {
    double slope, shift;
    slope = 2 / (upper_bound - lower_bound);
    shift = -1 * (upper_bound + lower_bound) / (upper_bound - lower_bound);
    return fmin(1, fmax(-1, slope * x + shift));
}