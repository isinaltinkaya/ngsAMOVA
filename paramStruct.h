#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__

#include "shared.h"

struct argStruct;

typedef struct paramStruct paramStruct;

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
    // TODO add argStruct ptr to paramStruct to avoid always passing both
    // argStruct *args;

    // number of sites non skipped for all individuals
    // nSites may not be !=totSites if minInd is set
    // or if a site is missing for all inds
    size_t nSites;
    // total number of sites processed
    size_t totSites;

    int ancder_nSites = 0;

    // ancestral allelic state for each site
    char *anc = NULL;
    // derived allelic state for each site
    char *der = NULL;

    int nInd;
    int nIndCmb;

    // TODO move this to amova analysis
    int nAmovaRuns;

    // input file type from enum
    int in_ft;
    char *DATETIME = NULL;

    // PRINT FUNCTIONS
    void printParams(FILE *fp);
    void printLut(FILE *fp);

    void init_LUTs();

    // validate that parameters make sense
    // e.g. nInd > 0
    void validate();

    void read_ancDerFile(char *fn);
};

/// @brief paramStruct_init initialize the paramStruct
/// @param args arguments argStruct
/// @return pointer to paramStruct
paramStruct *paramStruct_init(argStruct *args);
void paramStruct_destroy(paramStruct *p);

void check_consistency_args_pars(argStruct *args, paramStruct *pars);

#endif  // __PARAM_STRUCT__