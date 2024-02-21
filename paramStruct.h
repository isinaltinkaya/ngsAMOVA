#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__

#include "shared.h"

typedef struct strArray strArray;
typedef struct paramStruct paramStruct;
typedef struct formulaStruct formulaStruct;
typedef struct ibdStruct ibdStruct;

#define INPUT_IS_VCF \
    ( (pars->in_ft & IN_VCF) )

#define PROGRAM_NEEDS_METADATA \
    ( ( 0!= args->doAMOVA ) )

#define PROGRAM_NEEDS_FORMULA \
    ( ( 0!= args->doAMOVA ) )


/*
 * @typedef
 * @abstract paramStruct - parameter structure
 *
 * @field nSites				number of sites
 * @field nInd					number of individuals
 * @field pos					position
 *
 * @field nIndCmb				number of unique pairwise individual combinations
 *
 * @field major					major allele
 * @field minor					minor allele
 * @field ref					reference allele
 * @field anc					ancestral allele
 * @field der					derived allele
 */
struct paramStruct {

    // number of sites non skipped for all individuals
    // nSites may not be !=totSites if minInd is set
    // or if a site is missing for all inds
    size_t nSites;    // number of sites not-skipped for all individuals
    size_t totSites;  // total number of sites processed

    int nSites_arrays_size;
    int nContigs;
    int nInd;
    int nIndCmb;


    // input file type from enum
    int in_ft = 0;

    // a1: major/ref/ancestral allele
    // a2: minor/alt/derived allele
    // a1a2[site] = {base_a1,base_a2}
    // where base_* is 0 for A, 1 for C, 2 for G, 3 for T, 4 for BASE_UNOBSERVED
    int** a1a2;

    formulaStruct* formula;

    // ------------
    ibdStruct* ibd;
    // ------------

    int** pidx2inds;

    char* DATETIME;

    strArray* indNames;

};

/// @brief paramStruct_init initialize the paramStruct
/// @param args arguments argStruct
/// @return pointer to paramStruct
paramStruct* paramStruct_init(argStruct* args);
void paramStruct_destroy(paramStruct* p);

void check_consistency_args_pars(paramStruct* pars);


#endif  // __PARAM_STRUCT__
