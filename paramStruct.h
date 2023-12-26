#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__

#include "shared.h"

typedef struct paramStruct paramStruct;
typedef struct formulaStruct formulaStruct;
typedef struct ibdStruct ibdStruct;


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

    // a1: major/ref/ancestral allele
    // a2: minor/alt/derived allele
    // a1a2[site] = {base_a1,base_a2}
    // where base_* is 0 for A, 1 for C, 2 for G, 3 for T, 4 for BASE_UNOBSERVED
    int** a1a2 = NULL;

    int nSites_arrays_size;

    int nContigs;

    formulaStruct* formula = NULL;

    // ------------
    ibdStruct* ibd = NULL;
    // ------------

    int nInd;
    int nIndCmb;

    int** pidx2inds = NULL;

    // input file type from enum
    int in_ft = 0;

    char* DATETIME = NULL;

    // PRINT FUNCTIONS
    void printLut(FILE* fp);

    void init_LUTs();

    // validate that parameters make sense
    // e.g. nInd > 0
    void validate();

    char** indNames = NULL;

};

/// @brief paramStruct_init initialize the paramStruct
/// @param args arguments argStruct
/// @return pointer to paramStruct
paramStruct* paramStruct_init(argStruct* args);
void paramStruct_destroy(paramStruct* p);

void check_consistency_args_pars(paramStruct* pars);

// REQUIRE 
// in form require_object
// check if any of the specified analyses require object

bool require_formula(void);
bool require_metadata(void);
bool require_itemLabels(void);



#endif  // __PARAM_STRUCT__
