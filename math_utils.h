#ifndef __MATH_UTILS__
#define __MATH_UTILS__

#include "shared.h"



void gl_log10(int base, double errate, double *like);

void rescale_likelihood_ratio(double *like);
double log2ln(float ivar);


/*
 * Binomial coefficient
 * n choose k
 */
int nCk(int n, int k);

int nCk_idx(int nInd, int i1, int i2);

void prepare_LUT_indPair_idx(int nInd, int **LUT_indPair_idx);

namespace MATH {

	double SUM(double M[3][3]);
	double MEAN(double M[3][3]);
	double VAR(double M[3][3]);
	double SD(double M[3][3]);

	namespace EST{
		double Sij(double M[3][3]);
		double Fij(double M[3][3]);
	};

}	

#endif

