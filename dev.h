#include "vcfUtils.h"
#include "shared.h"
#include "em.h"




void DEV_input_VCF(argStruct *args, paramStruct *pars, formulaStruct *formulaSt, IO::outFilesStruct *outSt);

void* DEV_t_EM_2DSFS_GL3(void* p);

int DEV_EM_2DSFS_GL3(threadStruct* THREAD);

void DEV_prepare_distanceMatrix_originalData(argStruct *args, paramStruct *pars,  distanceMatrixStruct *dMS_orig, VCF::vcfData *VCF, pairStruct **pairSt,  formulaStruct *formulaSt, IO::outFilesStruct *outSt, sampleStruct *sampleSt);

void DEV_spawnThreads_pairEM_GL(argStruct *args, paramStruct *pars, pairStruct **pairSt, VCF::vcfData *VCF, IO::outFilesStruct *outSt, distanceMatrixStruct *distMatrix);
