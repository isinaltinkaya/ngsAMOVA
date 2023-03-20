#ifndef __EM__
#define __EM__

#include "shared.h"
#include "mathUtils.h"

#include <stdio.h>
#include <math.h>
#include <pthread.h>


void* t_EM_2DSFS_GL3(void* p);

int EM_2DSFS_GL3(threadStruct* THREAD);

void* DEV_t_EM_2DSFS_GL3(void* p);
int DEV_EM_2DSFS_GL3(threadStruct* THREAD);

void spawnThreads_pairEM(argStruct *args, paramStruct *pars, pairStruct **pairSt, vcfData *vcfd,  distanceMatrixStruct *distMatrix);

// void test_em(double **lngls3, double GLSFS[3][3],int **GTSFS,int i1,int i2,char *id1, char*id2, int pair_idx,size_t nSites, int shared_nSites, FILE *outfile);
// void test_em(double **lngls3, double *GLSFS,int **GTSFS,int i1,int i2,char *id1, char*id2, int pair_idx,size_t nSites, int shared_nSites, FILE *outfile);

#endif

