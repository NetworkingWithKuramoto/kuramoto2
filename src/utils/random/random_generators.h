#ifndef RANDOM_GENERATORS_H
#define RANDOM_GENERATORS_H

#include <stdint.h>

/**
 * @brief Generates uniform random numbers in the range [0, 1).
 * 
 * This function implements a portable random number generator (RNG) based on a lagged Fibonacci method.
 * 
 * @param[in,out] idum A pointer to an integer seed. The seed must be initialized with a negative value
 *                     on the first call, and it will be updated by the function on each subsequent call.
 * @return A uniform random float in the range [0, 1), exclusive of 1.
 */
float ran_lagged_fibonacci(int *idum);

/**
 * @brief Generates uniform random numbers in the range [0, 1).
 * 
 * This function implements a linear congruential generator (LCG) with additional shuffling for better randomness.
 * 
 * @param[in,out] idum A pointer to a long integer seed. The seed must be initialized with a negative value
 *                     on the first call, and it will be updated by the function on each subsequent call.
 * @return A uniform random float in the range [0, 1), exclusive of 1.
 */
float random_lin_cong(int *idum);
float generateRandomFloat(float min, float max);
double generate_gaussian(double mean, double variance);

float gasdev(int *idum);

/**
 * @brief Initializes the Mersenne Twister PRNG with a seed.
 * 
 * @param seed The initial seed for the generator.
 */
void init_mersenne_twister(uint32_t seed);

/**
 * @brief Generates a random unsigned 32-bit integer using Mersenne Twister.
 * 
 * @return A random unsigned 32-bit integer.
 */
uint32_t mersenne_twister_uint32(void);

/**
 * @brief Generates a random double in the range [0, 1) using Mersenne Twister.
 * 
 * @return A random double in the range [0, 1).
 */
double mersenne_twister_double(void);

/**
 * @brief Initializes the Xoshiro128+ PRNG with a seed.
 * 
 * @param seed The initial seed for the generator.
 */
void init_xoshiro128p(uint32_t seed);

/**
 * @brief Generates a random unsigned 32-bit integer using Xoshiro128+.
 * 
 * @return A random unsigned 32-bit integer.
 */
uint32_t xoshiro128p_uint32(void);

/**
 * @brief Generates a random double in the range [0, 1) using Xoshiro128+.
 * 
 * @return A random double in the range [0, 1).
 */
double xoshiro128p_double(void);

/**
 * @brief Initializes the PCG32 PRNG with a seed.
 * 
 * @param seed The initial seed for the generator.
 */
void init_pcg32(uint64_t seed);

/**
 * @brief Generates a random unsigned 32-bit integer using PCG32.
 * 
 * @return A random unsigned 32-bit integer.
 */
uint32_t pcg32_uint32(void);

/**
 * @brief Generates a random double in the range [0, 1) using PCG32.
 * 
 * @return A random double in the range [0, 1).
 */
double pcg32_double(void);

#endif // RANDOM_GENERATORS_H
