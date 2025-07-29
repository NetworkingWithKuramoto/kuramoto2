#include "random_generators.h"
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

#define IA 16807
#define IM 2147483647
#define AM (1.0 / (double)IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

float random_lin_cong(int *idum)
{
	/**
	* @brief Generates uniform random numbers in the range [0, 1).
	* 
	* This function implements a linear congruential generator (LCG) with additional shuffling for improved randomness.
	* It generates pseudo-random numbers based on a given seed `idum` and is designed for high statistical quality.
	* 
	* @param[in,out] idum A pointer to a long integer seed. This seed is updated on each function call to ensure randomness.
	*                     The seed must be initialized to a negative value before the first call.
	* 
	* @return A uniform random float in the range [0, 1), exclusive of 1.
	* 
	* @details
	* - The function uses the constants `IA`, `IM`, `IQ`, `IR`, `NDIV`, `NTAB`, and `EPS` to produce random numbers.
	* - The shuffling mechanism uses an array `iv` of size `NTAB` and an integer `iy` to further improve the randomness.
	* - If `idum` is initially non-positive or `iy` is zero, the array `iv` is initialized using the seed.
	* - The core of the random number generation is based on the linear congruential formula:
	*   ```
	*   idum = (IA * (idum % IQ)) - (IR * (idum / IQ));
	*   ```
	*   where `IA`, `IM`, `IQ`, and `IR` are constants chosen to achieve a long period and uniform distribution.
	* - The generated number is scaled to the range [0, 1) by multiplying by `AM` (which is `1.0/IM`).
	* - The function uses `RNMX` to ensure the random numbers do not reach exactly 1.
	* 
	* @note The function is reentrant and thread-safe as long as each thread uses its own seed.
	* @warning It is crucial to initialize `idum` to a negative value before the first function call.
	* 
	* ## Constants:
	* - `IA`, `IM`: Multiplier and modulus of the LCG.
	* - `IQ`, `IR`: Factors for the LCG division and subtraction.
	* - `NTAB`, `NDIV`: Size of the shuffle table and divisor for indexing into the table.
	* - `EPS`, `RNMX`: Machine precision and the maximum random number returned.
	*/

	int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }

    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0) *idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;

    // Ensure iy is within reasonable bounds
    iy = iy % IM;

    temp = AM * iy;

    if (temp > RNMX) return RNMX;
    else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software 7MZ9%"W5:!+). */

/**
* @brief Generates uniform random numbers in the range [0, 1).
* 
* This function implements a portable version of a random number generator (RNG) based on a lagged Fibonacci method. 
* It produces high-quality, uniform random numbers and is designed to minimize correlations over sequences.
* 
* @param[in,out] idum A pointer to an integer seed for initializing the random number generator.
*                     The first call should pass a negative value to properly initialize the sequence.
*                     Subsequent calls should use the positive, updated seed.
* 
* @return A uniform random float in the range [0, 1), exclusive of 1.
* 
* @details
* - The function uses a lagged Fibonacci generator with a period greater than `2^32`.
* - It initializes a shuffle table `ma` of size 56 using the seed `idum` and several constants:
*   - `MBIG`: The modulus used for the RNG.
*   - `MSEED`: A constant used to initialize the sequence.
*   - `MZ`: A zero value used for comparison.
*   - `FAC`: A scaling factor to map the output to the range [0, 1).
* - The `ma` array is filled and then scrambled based on several passes through its values.
* - After initialization, each call to `ran_lagged_fibonacci` updates two indices, `inext` and `inextp`, to produce a new random number.
*   The new number is calculated as the difference between two entries in the `ma` array.
* - The result is then scaled by `FAC` to produce a float in the desired range.
* 
* @note The function is reentrant and thread-safe as long as each thread uses its own seed.
* @warning To initialize the random number generator correctly, pass a negative value for `idum` on the first call.
* 
* ## Constants:
* - `MBIG`: The modulus used for generating the random numbers.
* - `MSEED`: The initial seed offset to avoid small seed effects.
* - `MZ`: The zero value to handle negative differences.
* - `FAC`: A scaling factor (`1.0 / MBIG`) to map random numbers to the range [0, 1).
*/
float ran_lagged_fibonacci(int *idum) {
    static int inext, inextp;
    static long ma[56];  // Store 55 values for the lagged Fibonacci generator
    static int initialized = 0;  // Track if initialization is done
    long mj, mk;
    int i, ii, k;

    // Initialize only on the first call OR when a negative seed is given
    if (*idum < 0 || !initialized) {
        initialized = 1;

        // Ensure positive seed for consistent behavior
        *idum = (*idum < 0) ? -(*idum) : *idum;

        // Initialize `mj` using modulus to prevent overflow issues
        mj = MSEED - (*idum % MBIG);
        if (mj < 0) mj += MBIG;
        ma[55] = mj;
        mk = 1;

        // Fill the `ma` array with shuffled values
        for (i = 1; i <= 54; i++) {
            ii = (21 * i) % 55;  // Shuffle pattern
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ) mk += MBIG;
            mj = ma[ii];
        }

        // Further randomize the `ma` table
        for (k = 1; k <= 4; k++) {
            for (i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        }

        inext = 0;
        inextp = 31;
    }

    // Generate the next random number
    if (++inext == 56) inext = 1;
    if (++inextp == 56) inextp = 1;

    mj = ma[inext] - ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext] = mj;

    return mj * FAC;  // Normalize result to [0,1)
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

void seedWithMicroseconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec); // Seed with microseconds
}

float generateRandomFloat(float min, float max) {
	seedWithMicroseconds();
    // Generate a random float between min and max
    float scale = rand() / (float) RAND_MAX; // Generates a float between 0 and 1
    return min + scale * (max - min);        // Scales it to the desired range
}

double generate_gaussian(double mean, double variance) {
	/**
	 * @brief Generates a random number with a Gaussian (normal) distribution.
	 *
	 * This function uses the Box-Muller transform to produce a random number
	 * with a Gaussian distribution, given a specified mean and variance.
	 * The Box-Muller transform converts two uniformly distributed random
	 * numbers into two independent standard normally distributed random
	 * numbers, which are then scaled to match the desired mean and variance.
	 *
	 * @param mean The mean (μ) of the Gaussian distribution.
	 * @param variance The variance (σ^2) of the Gaussian distribution.
	 * @return A random number following the Gaussian distribution
	 *         with the specified mean and variance.
	 *
	 * @note This function stores one of the generated Gaussian values for use
	 *       in the next function call, improving efficiency by avoiding
	 *       redundant computations.
	 */
    static int has_saved = 0;        // Flag to check if there's a saved number
    static double saved_value;       // Variable to store the saved number
    
    if (has_saved) {
        // If we have a saved number from a previous run, use it
        has_saved = 0;
        return mean + saved_value * sqrt(variance);
    }

    // Generate two uniform random numbers in the (0, 1) range
    double u1 = rand() / (double)RAND_MAX;
    double u2 = rand() / (double)RAND_MAX;

    // Apply Box-Muller Transform
    double radius = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;
    
    // Calculate the two Gaussian-distributed values
    double z0 = radius * cos(theta);
    double z1 = radius * sin(theta);
    
    // Save one for the next call
    saved_value = z1;
    has_saved = 1;

    // Transform the z0 value to match the desired mean and variance
    return mean + z0 * sqrt(variance);
}



float gasdev(int *idum)
{
	/**
	* @brief Generates normally distributed random numbers using the Box-Muller transform.
	* 
	* This function uses the Box-Muller transform to convert uniform random numbers into a normal distribution with mean 0 and standard deviation 1.
	* It employs a static variable `iset` to store one of the generated values, thus improving efficiency by returning two normal deviates per function call.
	* 
	* @param[in,out] idum A pointer to an integer seed for the uniform random number generator `ran_lagged_fibonacci`. 
	*                     The seed value is updated during the function's execution to maintain randomness.
	* 
	* @return A normally distributed random float with mean 0 and variance 1.
	* 
	* @details
	* - The function calls `ran_lagged_fibonacci` to generate uniform random numbers between -1 and 1.
	* - The generated numbers `v1` and `v2` are used to compute `rsq`, which is the sum of their squares.
	* - If `rsq` is greater than or equal to 1 or exactly 0, the process is repeated to find suitable `v1` and `v2`.
	* - Once valid `v1` and `v2` are obtained, the Box-Muller transformation is applied to produce two normal deviates.
	* - The result is stored in `gset` for the next function call, making the function efficient by using a static toggle `iset`.
	* - On subsequent calls, `gset` is returned directly, and `iset` is reset.
	* 
	* @note This function depends on `ran_lagged_fibonacci`, which is assumed to return uniform random numbers between 0 and 1.
	* @warning Ensure `idum` is initialized properly before calling this function to achieve randomness.
	*/

	float ran_lagged_fibonacci(int *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran_lagged_fibonacci(idum)-1.0;
			v2=2.0*ran_lagged_fibonacci(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

// ===========================================
//        Mersenne Twister Implementation
// ===========================================

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL
#define UPPER_MASK 0x80000000UL
#define LOWER_MASK 0x7fffffffUL

static uint32_t mt[N];
static int mti = N + 1;

void init_mersenne_twister(uint32_t seed) {
    mt[0] = seed;
    for (mti = 1; mti < N; mti++) {
        mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
    }
}

uint32_t mersenne_twister_uint32(void) {
    uint32_t y;
    static uint32_t mag01[2] = {0x0UL, MATRIX_A};

    if (mti >= N) {
        int kk;
        for (kk = 0; kk < N - M; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 1];
        }
        for (; kk < N - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 1];
        }
        y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 1];
        mti = 0;
    }

    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

double mersenne_twister_double(void) {
    return mersenne_twister_uint32() * (1.0 / 4294967296.0);
}

// ===========================================
//        Xoshiro128+ Implementation
// ===========================================

static uint32_t xoshiro_state[4];

void init_xoshiro128p(uint32_t seed) {
    xoshiro_state[0] = seed;
    xoshiro_state[1] = seed ^ 0x9e3779b9;
    xoshiro_state[2] = seed ^ 0x7f4a7c13;
    xoshiro_state[3] = seed ^ 0xf39c86b0;
}

uint32_t rotl(const uint32_t x, int k) {
    return (x << k) | (x >> (32 - k));
}

uint32_t xoshiro128p_uint32(void) {
    uint32_t result = xoshiro_state[0] + xoshiro_state[3];
    uint32_t t = xoshiro_state[1] << 9;

    xoshiro_state[2] ^= xoshiro_state[0];
    xoshiro_state[3] ^= xoshiro_state[1];
    xoshiro_state[1] ^= xoshiro_state[2];
    xoshiro_state[0] ^= xoshiro_state[3];

    xoshiro_state[2] ^= t;
    xoshiro_state[3] = rotl(xoshiro_state[3], 11);

    return result;
}

double xoshiro128p_double(void) {
    return xoshiro128p_uint32() * (1.0 / 4294967296.0);
}

// ===========================================
//        PCG32 Implementation
// ===========================================

static uint64_t pcg_state = 0x4d595df4d0f33173ULL;

void init_pcg32(uint64_t seed) {
    pcg_state = seed;
}

uint32_t pcg32_uint32(void) {
    pcg_state = pcg_state * 6364136223846793005ULL + 1;
    uint32_t xorshifted = ((pcg_state >> 18u) ^ pcg_state) >> 27u;
    uint32_t rot = pcg_state >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

double pcg32_double(void) {
    return pcg32_uint32() * (1.0 / 4294967296.0);
}
