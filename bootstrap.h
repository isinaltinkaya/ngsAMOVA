#ifndef __BOOTSTRAP__
#define __BOOTSTRAP__

#include "argStruct.h"
#include "dataStructs.h"
#include "shared.h"

/* FORWARD DECLARATIONS ----------------------------------------------------- */

typedef struct blockStruct blockStruct;
typedef struct blobStruct blobStruct;
typedef struct bootstrapData bootstrapData;
typedef struct bootstrapDataset bootstrapDataset;

/* -------------------------------------------------------------------------- */

/// @brief blockStruct - structure for storing a single block
/// @details
///   positions are 0-based
///   [start, end) - [inclusive start, exclusive end)
struct blockStruct {
    char *chr = NULL;
    int start = 0;  // inclusive
    int end = 0;    // exclusive
    int len = 0;    // end - start; length of the block
};

// to access the data
// if this block is not the last block in the contig
// 		loop: from this_blockStart to next_blockStart
// if this block is the last block in the contig
// 		loop: from this_blockStart to the end of the contig
// for (int gti = 0; gti < 3; gti++)
// {
// 	double x = vcfd->lngl[0][3 * vb + gti];
// }

/// @brief blobStruct - structure for storing all blocks
struct blobStruct {
    // Total number of blocks
    int nBlocks = 0;

    // /def blocks[nBlocks]
    // blocks[i] == pointer to a blockStruct at index i
    blockStruct **blocks = NULL;

    // /def blockPtrs[nBlocks]
    // blockPtrs[i] == pointer to the location of the data for block i
    // e.g. blockPtrs[42] == index of the first site in block 42 in vcf data
    int **blockPtrs = NULL;

    ~blobStruct();

    void addBlock();

    void _print();
};

// blobStruct *blobStruct_get(vcfData *vcf, argStruct *args);
blobStruct *blobStruct_get(vcfData *vcf, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt);
blobStruct *blobStruct_read_bed(const char *fn);
blobStruct *blobStruct_read_tab(const char *fn);
blobStruct *blobStruct_populate_blocks_withSize(vcfData *vcf, argStruct *args);

int sample_block_variant(metadataStruct *mtd, const int lvl, const int local_group_idx);

struct bootstrapData {
    // /def vblocks[nInd][nBlocks]
    // vblocks[i][j] == index of the vblock associated with the jth block of the ith synthetic individual
    //   == index of the individual in the original dataset that the jth block of the ith synthetic individual is copied from
    int **vblocks = NULL;

    // stats associated with the bootstrap
};

// bootstrapData *prepare_bootstrap_block_1level(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt, blobStruct *blobSt);
bootstrapData *prepare_bootstrap_block_1level(vcfData *vcfd, metadataStruct *mS, blobStruct *blobSt);
bootstrapData *prepare_bootstrap_blocks_multilevel(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt, blobStruct *blobSt);

struct bootstrapDataset {
    // /def bdata[nBootstraps]
    // bdata[i] == pointer to bootstrapData struct for ith bootstrap
    bootstrapData **bdata = NULL;

    int nBootstraps = 0;
    int nInd = 0;
    int nBlocks = 0;

    // stats for all bootstraps

    bootstrapDataset(int nBootstraps_, int nInd_, int nBlocks_);
    ~bootstrapDataset();
};

// bootstrapDataset* bootstrapDataset_get(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt,  blobStruct *blobSt);
bootstrapDataset *bootstrapDataset_get(vcfData *vcfd, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, blobStruct *blobSt);

#endif  // __BOOTSTRAP__