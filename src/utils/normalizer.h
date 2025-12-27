#ifndef NORMALIZER_H
#define NORMALIZER_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double norm_tanh(double x);
double norm_softsign(double x);
double norm_linear(double x, double a);
double norm_fcr(double x, double lower_bound, double upper_bound);
double wrap_to_pi(double x);

#endif