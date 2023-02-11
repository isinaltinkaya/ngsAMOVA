#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__

// #include "shared.h"

#include <stdio.h>

struct paramStruct;
struct argStruct;

/*
 * @typedef
 * @abstract paramStruct - parameter structure
 *
 * @field nSites				number of sites
 * @field nInd					number of individuals
 * @field pos					position
 *
 * @field lut_indsToIdx		lookup table for mapping two individuals
 * 								to their pair index
 * @field nIndCmb				number of unique pairwise individual combinations
 *
 * @field major					major allele
 * @field minor					minor allele
 * @field ref					reference allele
 * @field anc					ancestral allele
 * @field der					derived allele
 */
struct paramStruct
{

    // number of sites non skipped for all individuals
    // nSites may not be !=totSites if minInd is set
    // or if a site is missing for all inds
    size_t nSites;
    // total number of sites processed
    size_t totSites;

    int nInd;

    int **lut_indsToIdx;
    int **lut_idxToInds;

    int nIndCmb;

    int nAmovaRuns;
    int verbose;

    // input file type from enum
    int in_ft;

    char *DATETIME;

    // PRINT FUNCTIONS
    void printParams(FILE *fp);
    void printLut(FILE *fp);

    void init_LUTs();

    // validate that parameters make sense
    // e.g. nInd > 0
    void validate();

    // verbose print
    void vprint(const int verbose_threshold, const char *format, ...);
    void vprint(const char *format, ...);

};

/// @brief paramStruct_init initialize the paramStruct
/// @param args arguments argStruct
/// @return pointer to paramStruct
paramStruct *paramStruct_init(argStruct *args);
void paramStruct_destroy(paramStruct *p);

void check_consistency_args_pars(argStruct *args, paramStruct *pars);

#endif // __PARAM_STRUCT__