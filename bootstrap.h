#ifndef __BOOTSTRAP__
#define __BOOTSTRAP__

#include "amova.h"
#include "argStruct.h"
#include "dataStructs.h"
#include "paramStruct.h"
#include "shared.h"
#include "vcfReader.h"

/* FORWARD DECLARATIONS ----------------------------------------------------- */
struct distanceMatrixStruct;
struct metadataStruct;
struct vcfData;
struct amovaStruct;

typedef struct blockStruct blockStruct;
typedef struct blobStruct blobStruct;
typedef struct bootstrapReplicate bootstrapReplicate;
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

    void print();

    void get_bootstrap_replicates(vcfData *vcfd, bootstrapReplicate *bRep);
};

// blobStruct *blobStruct_get(vcfData *vcf);
blobStruct *blobStruct_get(vcfData *vcf, paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mS, formulaStruct *formulaSt);
blobStruct *blobStruct_get(vcfData *vcf, paramStruct *pars);
blobStruct *blobStruct_read_bed(const char *fn);
blobStruct *blobStruct_read_tab(const char *fn);
blobStruct *blobStruct_populate_blocks_withSize(vcfData *vcf);

void *blobStruct_destroy(blobStruct *blob);

int sample_block_variant(metadataStruct *mtd, const int lvl, const int local_group_idx);

bootstrapReplicate *get_bootstrap_replicate(vcfData *vcfd, blobStruct *blobSt);
int sample_with_replacement(int n);

struct bootstrapReplicate {
    // /def rBlocks[nBlocks] - set of block ids
    // sampled in the bootstrap replicate bootstrapReplicate
    int *rBlocks = NULL;

    distanceMatrixStruct *distanceMatrix = NULL;

    amovaStruct *amova = NULL;

    bootstrapReplicate(int nBlocks);
    ~bootstrapReplicate();
};

struct bootstrapDataset {
    // /def rep[nBootstraps]
    // rep[i] == pointer to bootstrapData struct for ith bootstrap replicate
    bootstrapReplicate **replicates = NULL;

    int nReplicates = 0;
    int nBlocks = 0;

    // \def phiValues[indexOfPhiStatistic][nBootstrapReplicates]
    double **phiValues = NULL;
    int nPhiValues = 0;

    void print_confidenceInterval(FILE *fp);

    void print();

    bootstrapDataset(paramStruct *pars, int nBootstraps_, int nBlocks_);
    ~bootstrapDataset();
};

bootstrapDataset *bootstrapDataset_get(vcfData *vcfd, paramStruct *pars, blobStruct *blobSt);

// distanceMatrixStruct *get_distanceMatrix_GL_bootstrapReplicate(const int nInd, const int nIndCmb, const int squareDistance, bootstrapReplicate *bRep, vcfData *vcfd, paramStruct *pars);
// distanceMatrixStruct *get_distanceMatrix_GT_bootstrapReplicate(const int nInd, const int nIndCmb, const int squareDistance, vcfData *vcfd, paramStruct *pars, blobStruct *blob, const int rep);
void get_distanceMatrix_GT(paramStruct *pars, distanceMatrixStruct *distanceMatrix, vcfData *vcfd, pairStruct **pairSt, blobStruct *blob);

#endif  // __BOOTSTRAP__