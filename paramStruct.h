#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__

#include "shared.h"

typedef struct allelesStruct allelesStruct;
typedef struct paramStruct paramStruct;
typedef struct formulaStruct formulaStruct;

typedef struct strArray{
	char** vals=NULL;
	int nvals;
	int nbuf;

	void add(const char* new_val);
	void print(FILE* fp);
	void print(void);
}strArray;
strArray* strArray_init(void);
void strArray_destroy(strArray* arr);



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

    int nContigs;

    allelesStruct *ancder = NULL;
    allelesStruct *majmin = NULL;
    formulaStruct *formula = NULL;

    int nInd;
    int nIndCmb;

    // input file type from enum
    int in_ft = 0;

    char *DATETIME = NULL;

    // PRINT FUNCTIONS
    void printLut(FILE *fp);

    void init_LUTs();

    // validate that parameters make sense
    // e.g. nInd > 0
    void validate();

    strArray* indNames = NULL;

};

/// @brief paramStruct_init initialize the paramStruct
/// @param args arguments argStruct
/// @return pointer to paramStruct
paramStruct *paramStruct_init(argStruct *args);
void paramStruct_destroy(paramStruct *p);

void check_consistency_args_pars(paramStruct *pars);

// REQUIRE 
// in form require_object
// check if any of the specified analyses require object

bool require_formula(void);
bool require_metadata(void);
bool require_itemLabels(void);



#endif  // __PARAM_STRUCT__
