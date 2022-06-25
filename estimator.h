#ifndef ESTIMATOR
#define ESTIMATOR

#include "shared.h"



void gl_log10(int base, double errate, double *like);

void rescale_likelihood_ratio(double *like);

double EM_2DSFS_GL10(double **lngl, double SFS[10][10], int i1, int i2, size_t nSites, double tole);
double EM_2DSFS_GL3(double **lngl, double SFS[3][3], int i1, int i2, size_t nSites, double tole,char *anc, char *der);

double log2ln(float ivar);

#endif

