#ifndef __ARG_STRUCT__
#define __ARG_STRUCT__


// #include "io.h"
#include "shared.h"

typedef struct argStruct argStruct;

extern u_char VERBOSE;
extern argStruct* args;

/*
 * @typedef
 * @abstract argStruct - argument structure
 *
 * @field *in_vcf_fn		input VCF filename
 * @field *in_dm_fn		input distance matrix filename
 * @field *in_mtd_fn	    input metadata filename
 * @field *in_blb_fn     input block bed filename
 * @field *in_dxy_fn     input dxy filename
 *
 * @field *out_fnp		pointer to output file prefix [angsdput]
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
 * @field doEM           [default=0]
 *                       [0] do NOT perform EM
 *                       [1] perform EM optimization with 3 GLs (MM,Mm,mm)
 *                               requires the ancestral and derived alleles to be specified
 *                       [2] perform EM optimization with 10 GLs
 *
 * @field doDxy          [default=0]
 *                       [0] do NOT perform Dxy
 *                       [1] perform Dxy for all pairs of groups in each hierarchical level
 *                               in the metadata file
 *                       ["str"] perform Dxy for all pairs of groups in hierarchical level "str" in the metadata file
 *                               e.g. "Regions" will perform Dxy for all pairs of groups in the "Regions"
 *                       ["str1,str2,str3"] perform Dxy for all pairs of groups in the given group list {str1,str2,str3}
 *                               e.g. "pop1,pop2,pop3" will perform Dxy for pop1-pop2, pop1-pop3, pop2-pop3.
 *
 * @field doDxyStr       [default=""] Stores the string given to doDxy, if any
 *                       [""] do NOT perform Dxy OR perform Dxy using numeric argument
 *
 * @field doPhylo           [default=0]
 *                       [0] do NOT perform neighbor joining
 *                       [1] perform neighbor joining using individuals as nodes
 *                       [2] perform neighbor joining using groups as nodes
 *
 * @field doDist		[0] use Sij similarity index
 * 						[1] use Dij (1-Sij) dissimilarity index
 *
 * @field sqDist		[0] use absolute value of distance measure (|dist_ij|)
 * 						[1] use squared distance measure (dist_ij^2)
 *
 * @field doInd			do ind pairs
 * @field ind1			ind1 id
 * @field ind2			ind2 id
 *
 *
 *
 *
 * @field nThreads		number of threads
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
 *                              [0] do NOT print distance matrix
 *                              [VALUE] print distance matrix
 *
 *
 * @field windowSize            [default = 0]
 *                              [0] do NOT use sliding window
 *                              [VALUE] use sliding window of size VALUE
 *
 *
 *
 */

struct argStruct {

    char* in_vcf_fn = NULL;
    char* in_dm_fn = NULL;
    char* in_mtd_fn = NULL;
    char* in_dxy_fn = NULL;

    char* out_fnp = NULL;

    char* in_region = NULL;
    char* in_regions_tab_fn = NULL;
    char* in_regions_bed_fn = NULL;

    int blockSize = 0;
    char* in_blocks_tab_fn = NULL;
    char* in_blocks_bed_fn = NULL;

    char* in_majorminor_fn = NULL;
    char* in_ancder_fn = NULL;

    char* command = NULL;

    char* formula = NULL;

    int doAMOVA = 0;
    int doEM = 0;
    int doDxy = 0;
    int doPhylo = 0;
    int doDist = 0;
    int doIbd = 0;

    int handle_neg_branch_length = 0;

    int printAmovaTable = 0;
    int printDistanceMatrix = 0;
    int printJointGenotypeCountMatrix = 0;
    int printBlocksTab = 0;

    int squareDistance = 1;

    char* dxyGroups = NULL;
    char* dxyLevels = NULL;

    // EM
    int minInd = -1;
    double tole = -1;
    int maxEmIter = -1;


    // ------------------------------------------------
    // VCF filters
    double min_af = 0.0;
    // ------------------------------------------------



    // ------------------------------------------------
    // IBDSEQ method implementation
    int ibdseq_minalleles = 2;
    double ibdseq_errormax = 0.001;
    double ibdseq_errorprop = 0.25;
    double ibdseq_ibdlod = 3.0;
    double ibdseq_ibdtrim = 0.0;



    int nThreads = 1;
    int seed = -1;
    int nBootstraps = 0;

    void check_arg_dependencies();

    // TODO check below
    char* doDxyStr = NULL;  // TODO maybe not needed
    int windowSize = 0;
    int* keyCols = NULL;
    int printDev = 0;

    void print(FILE* fp);
    void print();
};

//
/// @brief argStruct_get - parse command line arguments
/// @return pointer to the argStruct structure
argStruct* argStruct_get(int argc, char** argv);
void argStruct_destroy(argStruct* arg);

#endif  // __ARG_STRUCT__
