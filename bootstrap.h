#ifndef __BOOTSTRAP__
#define __BOOTSTRAP__



#include "shared.h"
#include "argStruct.h"
#include "dataStructs.h"

void prepare_bootstrap_blocks(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, sampleStruct *sampleSt, metadataStruct *mS, formulaStruct *formulaSt, IO::outFilesStruct *outSt, blobStruct *blobSt);


#endif // __BOOTSTRAP__