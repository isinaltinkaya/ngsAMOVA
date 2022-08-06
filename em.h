#ifndef __EM__
#define __EM__

#include "shared.h"
#include "math_utils.h"

#include <stdio.h>
#include <math.h>
#include <pthread.h>


void* t_EM_2DSFS_GL3(void* p);

int EM_2DSFS_GL3(IO::threadStruct* THREAD);

void test_em(double **lngls3, double GLSFS[3][3],int **GTSFS,int i1,int i2,char *id1, char*id2, int pair_idx,size_t nSites, int shared_nSites, FILE *outfile);

#endif

