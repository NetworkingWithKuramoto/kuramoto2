#include <math.h>
#include <stdlib.h>
#include "sampler.h"

void sample_logarithmically(int N, int M, double base, int start_from_zero, int samples[]) {
    for (int i = 0; i < M; i++) {
        double exponent = (double)i * log(N) / (log(base) * (M - 1));
        int value = (int)floor(pow(base, exponent));

        // Adjust for starting from zero or one
        if (start_from_zero) {
            value = (int)floor(pow(base, exponent) * (N - 1) / N); // Scale to go from 0 to N
        } else {
            value = (int)floor(pow(base, exponent));
        }

        if (value >= N) value = N; // Ensure last value is exactly N
        if (i > 0 && value <= samples[i - 1]) value = samples[i - 1] + 1; // Ensure uniqueness

        samples[i] = value;
    }

    // Force last element to be exactly N
    samples[M - 1] = N;
}