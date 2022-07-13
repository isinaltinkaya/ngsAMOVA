#ifndef __EM_OPTIMIZATION__
#define __EM_OPTIMIZATION__

#include "shared.h"
#include <stdio.h>


double EM_2DSFS_GL3(double **lngls, double SFS[3][3], int i1, int i2, size_t nSites, int shared_nSites, double tole, int *n_em_iter);

void test_em(double **lngls3, double GLSFS[3][3],int **GTSFS,int i1,int i2,char *id1, char*id2, int pair_idx,size_t nSites, int shared_nSites, FILE *outfile);

#endif

