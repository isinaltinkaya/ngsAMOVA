#include "vcfUtils.h"
#include "shared.h"
#include "em.h"

void DEV_input_VCF(argStruct *args, paramStruct *pars, formulaStruct *formulaSt, IO::outFilesStruct *outSt);

void *DEV_t_EM_2DSFS_GL3(void *p);

int DEV_EM_2DSFS_GL3(threadStruct *THREAD);

void DEV_prepare_distanceMatrix_originalData(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS_orig, vcfData *vcfd, pairStruct **pairSt, formulaStruct *formulaSt, IO::outFilesStruct *outSt, sampleStruct *sampleSt);

void DEV_spawnThreads_pairEM_GL(argStruct *args, paramStruct *pars, pairStruct **pairSt, vcfData *vcfd, IO::outFilesStruct *outSt, distanceMatrixStruct *distMatrix);
