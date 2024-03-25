#ifndef __MATH_UTILS__
#define __MATH_UTILS__

#include <math.h>
#include <stdio.h>

#include <limits>

#include "shared.h"
#include "vcfReader.h"

#define SQUARE(n) ((((n)) * ((n))))

/*
 * Macro:[LOG2LN]
 * convert log_10(x) to log_e(x)
 */
#define LOG2LN(x) ((((x)) / (M_LOG10E)))

#define LN2LOG(x) ((((x) * (M_LOG10E))))

/*
 * Binomial coefficient: n choose k
 *
 * @brief nCk - n choose k recursive function
 * @param n
 * @param k
 * @return choose(n,k)
 */

inline int nCk(int n, int k) {
    if (k == 0) {
        return 1;
    }
    return (n * nCk(n - 1, k - 1)) / k;
}


/*
 * [nCk_idx]
 * Maps a given pair of objects to their index in
 * the lexicographically ordered binomial coefficients
 * a.k.a. array of pairs
 *
 * @param nInd: number of individuals
 * @param i1: index of first individual
 * @param i2: index of second individual
 *
 * @return: index of pair in lexicographically ordered array of pairs
 *
 * @example:
 *
 *  	0	1	2	3
 * 0		01	02	03
 * 1			12	13
 * 2				23
 * 3
 *
 * 01 = 0, 02 = 1, 03 = 2, 12 = 3, 13 = 4, 23 = 5
 *
 * Formula:
 * (nCk(nInd, 2) - nCk((nInd - i1), 2)) + (i2 - i1) - 1;
 *
 * e.g.:
 * nInd=5, i1=1, i2=3
 * nCk(5,2) - nCk(4,2) + (3-1) - 1 = 5 - 6 + 2 - 1 = 0
 *
 * Total number of pairs: nCk(nInd, 2)
 * Starting index for each individual: nCk(nInd, 2) - nCk((nInd - i), 2)
 * Difference between i2 and i1: (i2 - i1)
 * Index of pair: (nCk(nInd, 2) - nCk((nInd - i1), 2)) + (i2 - i1) - 1;
 *
 */

inline int nCk_idx(int nInd, int i1, int i2) {
    ASSERT(i1 < nInd && i2 < nInd);  // safeguard for wrong order of arguments
    if (i2 > i1) {
        return (nCk(nInd, 2) - nCk((nInd - i1), 2)) + (i2 - i1) - 1;
    } else {
        return (nCk(nInd, 2) - nCk((nInd - i2), 2)) + (i1 - i2) - 1;
    }
}


/*
 * rand()
 * 		generates a random integer between 0 (incl) and RAND_MAX (incl)
 *
 *
 * rand() / (RAND_MAX + 1.0)
 * 		generates a floating point number [0,1)
 * 		max value is RAND_MAX / (RAND_MAX + 1.0)
 * 		which is 0.99996... (for RAND_MAX = 32767)
 * 		guaranteed min value for RAND_MAX is 32767
 *
 */



#endif  // __MATH_UTILS__
