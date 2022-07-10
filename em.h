#ifndef __EM_OPTIMIZATION__
#define __EM_OPTIMIZATION__

#include "shared.h"

// int EM_2DSFS_GL3(double **lngl, double SFS[3][3], int i1, int i2, size_t nSites, int shared_nSites,double tole,char *anc, char *der);
// int EM_2DSFS_GL3(double **lngl, double SFS[3][3], int i1, int i2, size_t nSites, int shared_nSites, double tole);

double EM_2DSFS_GL3(double **lngls, double SFS[3][3], int i1, int i2, size_t nSites, int shared_nSites, double tole, int *n_em_iter);

#endif

