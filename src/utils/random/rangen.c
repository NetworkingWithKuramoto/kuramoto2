#include "rangen.h"

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran31(int *idum)
{
	/**
	* @brief Generates uniform random numbers in the range [0, 1).
	* 
	* These functions (`ran31` and `ran32`) implement a random number generator (RNG) using the lagged Fibonacci method.
	* They are designed to produce high-quality uniform random numbers over a sequence and are based on the same logic.
	* The difference between `ran31` and `ran32` may lie in their intended use in different contexts of the program, 
	* but both generate random numbers using the same algorithm.
	* 
	* @param[in,out] idum A pointer to an integer seed for initializing the random number generator.
	*                     The first call should pass a negative value to properly initialize the sequence.
	*                     Subsequent calls should use the positive, updated seed.
	* 
	* @return A uniform random double in the range [0, 1), exclusive of 1.
	* 
	* @details
	* - The functions use a lagged Fibonacci generator to produce random numbers.
	* - The internal array `ma` of size 56 is initialized and scrambled based on the seed `idum`.
	*   - `MBIG`: The modulus used for the RNG.
	*   - `MSEED`: A large constant to avoid small seed issues.
	*   - `MZ`: A zero value to handle negative differences.
	*   - `FAC`: A scaling factor to map the output to the range [0, 1).
	* - On initialization (when `*idum` is negative or the first call), the shuffle table `ma` is filled and scrambled to ensure randomness.
	* - Each call to the functions updates indices `inext` and `inextp` to produce a new random number from the table.
	* - The generated random number is scaled by `FAC` to produce a double in the range [0, 1).
	* 
	* @note
	* - The functions are reentrant and thread-safe as long as each thread uses its own seed.
	* - To correctly initialize the random number generator, pass a negative value for `idum` on the first call.
	* 
	* ## Constants
	* - `MBIG`: The modulus value used for generating the random numbers.
	* - `MSEED`: The initial seed offset to avoid small seed effects.
	* - `MZ`: The zero value for ensuring positive random values.
	* - `FAC`: A scaling factor (`1.0 / MBIG`) to map the output random numbers to the range [0, 1).
	*/

	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran32(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

