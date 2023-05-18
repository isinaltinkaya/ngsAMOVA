#ifndef __EM__
#define __EM__

#include <pthread.h>

#include "argStruct.h"
#include "mathUtils.h"

void* DEV_t_EM_2DSFS_GL3(void* p);
int DEV_EM_2DSFS_GL3(indPairThreads* THREAD);

/// @brief spawnThreads_pairEM spawn threads for running EM algorithm for each individual pair
void spawnThreads_pairEM(paramStruct* pars, pairStruct** pairSt, vcfData* vcfd, distanceMatrixStruct* distMatrix);

/// @brief thread handler for EM_optim_jgd_gl3
void* t_EM_optim_jgd_gl3(void* p);

/// @brief EM algorithm for 3x3 pairwise genotype categories
int EM_optim_jgd_gl3(indPairThreads* THREAD);

#endif  // __EM__