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

#endif

