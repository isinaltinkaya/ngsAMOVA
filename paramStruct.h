#ifndef __PARAM_STRUCT__
#define __PARAM_STRUCT__

#include "shared.h"

struct argStruct;

typedef struct alleleStruct alleleStruct;

typedef struct paramStruct paramStruct;

typedef struct formulaStruct formulaStruct;

struct formulaStruct {
    // @nTokens number of tokens in the formula
    // e.g. formula: "Individual ~ Region/Population/Subpopulation"
    // 		nTokens = 4
    // 		corresponds to 3 hierarchical levels (Region, Population, Subpopulation)
    // 		thus nTokens == nLevels + 1
    int nTokens = 0;

    // the formula in the raw text form as it is in the argument
    // e.g. "Individual ~ Region/Population/Subpopulation"
    char *formula = NULL;

    // @formulaTokens
    // array of tokens in the formula
    // e.g. formula: "Individual ~ Region/Population/Subpopulation"
    // 		formulaTokens = {"Individual","Region","Population","Subpopulation"}
    char **formulaTokens;

    // @formulaTokenIdx[nTokens]
    //
    // maps index of the token in formula to index of the corresponding column in metadata file
    // formulaTokenIdx[indexOfTokenInFormula] = indexOfTokenInMetadataFile
    //
    // e.g. metadata file header: "Individual,Population,Etc,Region,Subpopulation"
    // 		formula: "Individual ~ Region/Population/Subpopulation"
    // 		formulaTokenIdx = {0,3,1,2}
    int *formulaTokenIdx;

    void print(FILE *fp);

    /// match the given metadata token with formula tokens
    /// @param mtd_tok 		- metadata token to match
    /// @param mtd_col_idx	- index of the metadata column containing mtd_tok
    /// @return int			- index if found any match, -1 otherwise
    int setFormulaTokenIdx(const char *mtd_tok, const int mtd_col_idx);

    // TODO deprec
    //  @brief shrink - shrink the size of the arrays defined with default max values to the actual size needed
    void shrink();
};

formulaStruct *formulaStruct_get(const char *formula);
void formulaStruct_validate(formulaStruct *fos, const int nLevels);
void formulaStruct_destroy(formulaStruct *fos);
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
    size_t totSites;  // total number of sites processed

    int nContigs;

    alleleStruct *ancder = NULL;
    alleleStruct *majmin = NULL;
    formulaStruct *formula = NULL;

    int nInd;
    int nIndCmb;

    // TODO move this to amova analysis
    int nAmovaRuns = 0;

    // input file type from enum
    int in_ft = 0;
    char *DATETIME = NULL;

    // PRINT FUNCTIONS
    void printParams(FILE *fp);
    void printLut(FILE *fp);

    void init_LUTs();

    // validate that parameters make sense
    // e.g. nInd > 0
    void validate();
};

struct alleleStruct {
    // \def nSites[nContigs]
    //      nSites[contig_i] == number of sites with a1 a2 data in contig_i
    int *nSites = NULL;

    int nContigs = 0;

    // \def a1[nContigs][nSites[contig_i]]
    //      a1[contig_i][site_j] == a1 data for site_j in contig_i
    //      a1 can be ancestral or major allele
    char **a1 = NULL;

    // \def a2[nContigs][nSites[contig_i]]
    //      a2[contig_i][site_j] == a2 data for site_j in contig_i
    //      a2 can be derived or minor allele
    char **a2 = NULL;

    alleleStruct();
    ~alleleStruct();
};

/// @brief alleleStruct_read read alleles file
/// @param fn alleles file name
/// @return pointer to alleleStruct
/// @details alleles files are tab delimited files with 4 columns:
/// 1. contig name
/// 2. position
/// 3. a1
/// 4. a2
/// a1 and a2 can be ancestral and derived or major and minor alleles
alleleStruct *alleleStruct_read(const char *fn);

/// @brief alleleStruct_destroy destroy alleleStruct
/// @param as pointer to alleleStruct
void alleleStruct_destroy(alleleStruct *A);

/// @brief paramStruct_init initialize the paramStruct
/// @param args arguments argStruct
/// @return pointer to paramStruct
paramStruct *paramStruct_init(argStruct *args);
void paramStruct_destroy(paramStruct *p);

void check_consistency_args_pars(argStruct *args, paramStruct *pars);

#endif  // __PARAM_STRUCT__