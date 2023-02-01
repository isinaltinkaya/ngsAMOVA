#ifndef __mathUtils__
#define __mathUtils__

#include "shared.h"
#include "vcfUtils.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits>




#define SQUARE(n) ((n)*(n))



/*
 * Macro:[LOG2LN]
 * convert log_10(x) to log_e(x)
 */
#define LOG2LN(x) (x/M_LOG10E)


#define LN2LOG(x) (x*M_LOG10E)


void gl_log10(int base, double errate, double *like);



/*
 * Binomial coefficient
 * n choose k
 */
int nCk(int n, int k);


/*
 * Maps a given pair of objects to their index in
 * the lexicographically ordered binomial coefficients
 * a.k.a. array of pairs
 */
int nCk_idx(int nInd, int i1, int i2);


/* [Look-Up Table]
 * Prepare individual pair - pair index lookup tables
 *
 */
void set_lut_indsToIdx_2way(int nInd, int nIndCmb, int** lut_indsToIdx, int** lut_idxToInds);

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

// int sample_

namespace MATH {

	double SUM(double M[3][3]);
	double SUM(double* M);
	double MEAN(double M[3][3]);
	double MEAN(double* M);
	double VAR(double M[3][3]);
	double VAR(double* M);
	double SD(double M[3][3]);
	double SD(double* M);
	
	double SUM(int* M);
	double MEAN(int* M);
	double VAR(int* M);
	double SD(int* M);


	//ESTIMATE?
	namespace EST{

		double Sij(double M[3][3]);
		double Fij(double M[3][3]);
		double IBS0(double M[3][3]);
		double IBS1(double M[3][3]);
		double IBS2(double M[3][3]);
		double R0(double M[3][3]);
		double R1(double M[3][3]);
		double Kin(double M[3][3]);
		
		double Sij(double* M);
		double Fij(double* M);
		double IBS0(double* M);
		double IBS1(double* M);
		double IBS2(double* M);
		double R0(double* M);
		double R1(double* M);
		double Kin(double* M);

		double Sij(int* M,int S);
		double Fij(int* M, int S);
		double IBS0(int* M, int S);
		double IBS1(int* M, int S);
		double IBS2(int* M, int S);
		double R0(int* M, int S);
		double R1(int* M, int S);
		double Kin(int* M, int S);


		double Dij(int *M, int S);
		double Dij(double *M);
		double Dij(double M[3][3]);
	};

}	


#endif

