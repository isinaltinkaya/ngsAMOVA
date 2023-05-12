#ifndef __MATH_UTILS__
#define __MATH_UTILS__

#include <math.h>
#include <stdio.h>

#include <limits>

#include "shared.h"
#include "vcfReader.h"

#define SQUARE(n) ((n) * (n))

/*
 * Macro:[LOG2LN]
 * convert log_10(x) to log_e(x)
 */
#define LOG2LN(x) (x / M_LOG10E)

#define LN2LOG(x) (x * M_LOG10E)

void gl_log10(int base, double errate, double *like);

/*
 * Binomial coefficient: n choose k
 *
 * @brief nCk - n choose k recursive function
 * @param n
 * @param k
 * @return choose(n,k)
 */
int nCk(int n, int k);

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
int nCk_idx(int nInd, int i1, int i2);

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

namespace MATH {

double SUM(double M[3][3]);
double SUM(double *M);
double MEAN(double M[3][3]);
double MEAN(double *M);
double VAR(double M[3][3]);
double VAR(double *M);
double SD(double M[3][3]);
double SD(double *M);

double SUM(int *M);
double MEAN(int *M);
double VAR(int *M);
double SD(int *M);

/*
 *
 * [S_ij Similarity index]
 *
 * Similarity index (S_ij) based on the
 * probability of identity by descent
 *
 * 00 01 02 10 11 12 20 21 22
 * A  D  G  B  E  H  C  F  I
 *
 * A + I + ((B+D+E+F+H)/2)
 *
 */

/*
 * Reference for Fij IBS0 IBS1 IBS2 R0 R1 Kin
 * https://doi.org/10.1111/mec.14954
 */

/*
 *
 * [F_ij F statistic]
 *
 * 00 01 02 10 11 12 20 21 22
 * A  D  G  B  E  H  C  F  I
 *
 * (2C+2G-E) / (2C+2G+B+D+E+F+H)
 *
 */
double Sij(double M[3][3]);
double Sij(double *M);
double Sij(int *M, int S);

double Fij(double *M);
double Fij(int *M, int S);
double Fij(double M[3][3]);

double IBS0(double *M);
double IBS0(int *M, int S);
double IBS0(double M[3][3]);

double IBS1(double *M);
double IBS1(int *M, int S);
double IBS1(double M[3][3]);

double IBS2(double *M);
double IBS2(int *M, int S);
double IBS2(double M[3][3]);

double R0(double *M);
double R0(int *M, int S);
double R0(double M[3][3]);

double R1(double *M);
double R1(int *M, int S);
double R1(double M[3][3]);

double Kin(double *M);
double Kin(int *M, int S);
double Kin(double M[3][3]);

double Dij(int *M, int S);
double Dij(double *M);
double Dij(double M[3][3]);

}  // namespace MATH

#endif  // __MATH_UTILS__
