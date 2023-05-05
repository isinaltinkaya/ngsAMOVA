#ifndef __EM__
#define __EM__

#include <pthread.h>

#include "mathUtils.h"

void* t_EM_optim_jgd_gl3(void* p);

int EM_optim_jgd_gl3(threadStruct* THREAD);

void* DEV_t_EM_2DSFS_GL3(void* p);
int DEV_EM_2DSFS_GL3(threadStruct* THREAD);

void spawnThreads_pairEM(argStruct* args, paramStruct* pars, pairStruct** pairSt, vcfData* vcfd, distanceMatrixStruct* distMatrix);

#endif  // __EM__