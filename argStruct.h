#ifndef __ARG_STRUCT__
#define __ARG_STRUCT__

#include "shared.h"
#include "ctype.h"

// TODO is this necessary?
struct argStruct;

/*
 * @typedef
 * @abstract argStruct - argument structure
 *
 * @field *in_vcf_fn		pointer to input file name
 * @field *in_mtd_fn	pointer to input Metadata file name
 * @field *out_fn		pointer to output file prefix [angsdput]
 *
 * @field isSim			input is vcfgl simulation output
 * 							anc=ref and der=alt[0]
 *
 * @field isTest		[DEV] run for testing purposes
 *
 * @field minInd		[-1 = not set]
 * 						minimum number of individuals needed
 * 						for site to be included in analyses
 *
 * 						if minInd not set; set minInd=2
 * 						== no filter, include all sites that exist in pair
 *
 * 						if minInd==nInd; set minInd=0
 * 						== site should be nonmissing for all individuals
 * 						to be included
 *
 * @field blockSize		[0 = not set]
 *						Block size to be used in block bootstrapping
 *
 *
 * @field formula		[not set]
 * 						the formula of the AMOVA model to be fitted
 *
 * 						of the form:
 *  						{LEFT_HAND_SIDE} ~ {RIGHT_HAND_SIDE}
 *
 *	 						LEFT_HAND_SIDE: Column defining the distance matrix
 *											i.e. Samples
 *
 *	 						RIGHT_HAND_SIDE: Column(s) defining the hierarchical stratification levels
 *											i.e. Populations, Regions, etc
 *
 *
 * 						e.g. 1 level
 * 							Samples ~ Populations
 * 							{SamplesColumn} ~ {HierLvl1_StrataColumn}
 *
 * 						e.g. 2 levels
 * 							Samples ~ Regions/Populations
 * 							{SamplesColumn} ~ {HierLvl1_StrataColumn} / {HierLvl2_StrataColumn}
 *
 * 						e.g. 3 levels
 * 							Samples ~ Continents/Regions/Populations
 * 							{SamplesColumn} ~ {HierLvl1_StrataColumn} / {HierLvl2_StrataColumn} / {HierLvl3_StrataColumn}
 *
 * 						NOTE: All column names must exist in Metadata file
 * 						e.g.
 *
 * 						Samples,Populations,Regions
 * 						ind1,Pop1,Reg1
 * 						ind2,Pop1,Reg1
 * 						ind3,Pop2,Reg2
 * 						ind4,Pop2,Reg2
 *
 *
 * @field keyCols		defines the index(es) (1-based) of the
 * 						column(s) in Metadata file 'in_mtd_fn'
 * 						to be used to define hierarchical stratification levels
 *
 * @field doAMOVA		[0]
 * 						1 use 10 genotype likelihoods (GL)
 * 						2 use genotypes (GT) (NOTE: Only for benchmark purposes)
 *
 * @field doDist		[0] use Sij similarity index
 * 						[1] use Dij (1-Sij) dissimilarity index
 * 						[2] use Fij F statistic [DEPRECATED]
 *
 * @field sqDist		[0] use absolute value of distance measure (|dist_ij|)
 * 						[1] use squared distance measure (dist_ij^2)
 *
 * @field doInd			do ind pairs
 * @field ind1			ind1 id
 * @field ind2			ind2 id
 *
 *
 * @field doTest		test for convergence
 *
 *
 * @field mThreads		maximum number of threads defined by user
 * @field mEmIter		maximum number of iterations allowed for em
 *
 * @field seed			random seed for bootstrapping
 *
 * @field nBootstraps	number of bootstraps to be performed for AMOVA significance test
 *
 * 						bootstraps[0] = original data
 * 						bootstraps[1..nBootstraps] = bootstrapped data
 * 							size is nBootstraps+1
 *
 *
 * @field printMatrix           [default = 0]
 *                      [0] do NOT print distance matrix
 *                      [1] print distance matrix in human-readable format
 *                      [2] print distance matrix in gzipped format
 *
 * @field windowSize            [default = 0]
 *                              [0] do NOT use sliding window
 *                              [VALUE] use sliding window of size VALUE
 * 
 * @field printJointGenoDist    [default = 0]
 *                              [0] do NOT print joint genotype distributions of pairs of individuals
 */
struct argStruct
{

    int verbose;

    char *in_vcf_fn;
    char *in_dm_fn;
    char *in_mtd_fn;
    char *in_jgcd_fn;
    char *in_jgpd_fn;
    char *out_fn;

    char *command;

    char *formula;
    int *keyCols;

    int blockSize;
    int doAMOVA;
    int doEM;

    int printAmovaTable;
    int printMatrix;
    int printJointGenoCountDist;
    int printJointGenoProbDist;
    int printDev;

    int isSim;
    int isTest;
    int doDist;
    int do_square_distance;
    int minInd;


    int hasColNames;

    double tole;
    int doTest;

    int mThreads;
    int mEmIter;

    int seed;
    int nBootstraps;

    int gl2gt;

    
    int windowSize;
};

//
/// @brief argStruct_init - initialize the argStruct structure
/// @return pointer to the argStruct structure
argStruct *argStruct_init();

argStruct *argStruct_get(int argc, char **argv);
void argStruct_destroy(argStruct *arg);
void argStruct_print(FILE *fp, argStruct *arg);

#endif // __ARG_STRUCT__