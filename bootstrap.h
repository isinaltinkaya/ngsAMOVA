#ifndef __BOOTSTRAP__
#define __BOOTSTRAP__



#include "shared.h"
#include "argStruct.h"
#include "dataStructs.h"


int sample_block_variant(metadataStruct* mtd, const int r_lvl_idx);

typedef struct bootstrapData{

	// /def vblocks[nInd][nBlocks]
	// vblocks[i][j] == index of the vblock associated with the jth block of the ith synthetic individual
	//   == index of the individual in the original dataset that the jth block of the ith synthetic individual is copied from
    int** vblocks = NULL;

    // stats associated with the bootstrap


} bootstrapData;

bootstrapData *prepare_bootstrap_block_1level(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt, blobStruct *blobSt);
bootstrapData *prepare_bootstrap_blocks_multilevel(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt, blobStruct *blobSt);

typedef struct bootstrapDataset{

	// /def bdata[nBootstraps]
	// bdata[i] == pointer to bootstrapData struct for ith bootstrap
    bootstrapData **bdata = NULL;

    int nBootstraps=0;
    int nInd=0;
    int nBlocks=0;

    // stats for all bootstraps

    bootstrapDataset(int nBootstraps_, int nInd_, int nBlocks_);
    ~bootstrapDataset();

} bootstrapDataset;

bootstrapDataset* bootstrapDataset_get(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt,  blobStruct *blobSt);


#endif // __BOOTSTRAP__