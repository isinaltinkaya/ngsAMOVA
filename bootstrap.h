#ifndef __BOOTSTRAP__
#define __BOOTSTRAP__

#include "argStruct.h"
#include "dataStructs.h"
#include "shared.h"
// #include "vcfReader.h"

/* FORWARD DECLARATIONS ----------------------------------------------------- */

typedef struct blockStruct blockStruct;
typedef struct blobStruct blobStruct;
typedef struct bootstrapRep bootstrapRep;
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

    bootstrapDataset *bootstraps = NULL;

    ~blobStruct();

    void addBlock();

    void _print();

    void print(IO::outputStruct *out_dm_fs);

    void get_bootstrap_replicates(vcfData *vcfd, bootstrapRep *bRep);
};

// blobStruct *blobStruct_get(vcfData *vcf, argStruct *args);
blobStruct *blobStruct_get(vcfData *vcf, paramStruct *pars, argStruct *args, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt);
blobStruct *blobStruct_get(vcfData *vcf, paramStruct *pars, argStruct *args);
blobStruct *blobStruct_read_bed(const char *fn);
blobStruct *blobStruct_read_tab(const char *fn);
blobStruct *blobStruct_populate_blocks_withSize(vcfData *vcf, argStruct *args);

void *blobStruct_destroy(blobStruct *blob);

int sample_block_variant(metadataStruct *mtd, const int lvl, const int local_group_idx);

bootstrapRep *get_bootstrap_replicate(vcfData *vcfd, blobStruct *blobSt);
int sample_with_replacement(int n);

struct bootstrapRep {
    // /def rBlocks[nBlocks] - set of block ids
    // sampled in the bootstrap replicate bootstrapRep
    int *rBlocks = NULL;

    // stats associated with the bootstrap

    bootstrapRep(int nBlocks);
    ~bootstrapRep();
};

struct bootstrapDataset {
    // /def bootRep[nBootstraps]
    // bootRep[i] == pointer to bootstrapData struct for ith bootstrap replicate
    bootstrapRep **bootRep = NULL;

    // /def distanceMatrixRep[nBootstraps]
    // distanceMatrixRep[i] == pointer to distanceMatrixStruct for ith bootstrap replicate
    distanceMatrixStruct **distanceMatrixRep = NULL;

    int nBootstraps = 0;
    int nBlocks = 0;

    // stats for all bootstraps

    bootstrapDataset(argStruct *args, paramStruct *pars, int nBootstraps_, int nBlocks_);
    ~bootstrapDataset();
};

bootstrapDataset *bootstrapDataset_get(vcfData *vcfd, paramStruct *pars, argStruct *args, blobStruct *blobSt);
distanceMatrixStruct *get_distanceMatrix_GL_bootstrapRep(const int nInd, const int nIndCmb, const int squareDistance, bootstrapRep *bRep, vcfData *vcfd);
distanceMatrixStruct *get_distanceMatrix_GT_bootstrapRep(const int nInd, const int nIndCmb, const int squareDistance, bootstrapRep *bRep, vcfData *vcfd);
void readSites_bootstrapReplicate(vcfData *vcfd, argStruct *args, paramStruct *pars, pairStruct **pairSt, bootstrapRep *brep);

#endif  // __BOOTSTRAP__