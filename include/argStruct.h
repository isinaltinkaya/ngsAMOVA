#ifndef __ARG_STRUCT__
#define __ARG_STRUCT__

#include "dmat.h"


typedef struct argStruct argStruct;

extern uint8_t PROGRAM_VERBOSITY_LEVEL;
extern char* PROGRAM_VERSION_INFO;
extern char* PROGRAM_COMMAND;
extern argStruct* args;

struct argStruct {

    char* in_vcf_fn = NULL;
    char* in_dm_fn = NULL;
    char* in_mtd_fn = NULL;
    char* in_jgtmat_fn = NULL;

    char* out_fnp = NULL;

    char* in_region = NULL;
    char* in_regions_tab_fn = NULL;
    char* in_regions_bed_fn = NULL;

    int blockSize = -1;
    char* in_blocks_tab_fn = NULL;
    char* in_blocks_bed_fn = NULL;

    char* in_majorminor_fn = NULL;

    char* in_a2f_fn = NULL;

    char* command = NULL;

    // the formula in the raw text form as it is in the argument
    // e.g. "Individual ~ Region/Population/Subpopulation"
    char* formula = NULL;


    int doJGTM = 0;
    int doDist = 0;
    int doEM = 0;
    int doAMOVA = 0;
    int doDxy = 0;
    int doPhylo = 0;
    int doIbd = 0;
    int doMajorMinor = 0;
    int doBlockBootstrap = 0;

    int doUnitTests = 0;

    //TODO implement
    int doDryRun = 0;

    int bcfSrc = 0;

    int amova_euclid = 1;


    int alloc_strategy = ARG_ALLOC_STRATEGY_BASIC;


    int print_jgtm = -1;
    int print_jgtm_ctype = -1;
    int print_dm = -1;
    int print_dm_ctype = -1;
    int print_pruned_dm = -1;
    int print_pruned_dm_ctype = -1;
    int print_amova = -1;
    int print_amova_ctype = -1;
    int print_blocks = -1;
    int print_blocks_ctype = -1;
    int print_bootstrap = -1;
    int print_bootstrap_ctype = -1;
    int print_tree = -1;
    int print_tree_ctype = -1;
    int print_dxy = -1;
    int print_dxy_ctype = -1;
    int print_ibd = -1;
    int print_ibd_ctype = -1;



    // ------------------------------------------------
    // block bootstrapping
    int nBootstraps = 0;
    double bootstrap_pctci = -1.0;

    // doPhylo
    int handle_neg_branch_length = 0;


    // EM
    double tole = -1;
    int maxEmIter = -1;

    // dmat
    int dm_method = -1;
    int dm_transform = -1;

    int prune_dmat = 0;


    // ------------------------------------------------
    // VCF filters
    double min_info_a2f = 0.0;
    int minInd = -1;
    int min_a2c = 0;
    // ------------------------------------------------


    double* a2freqs = NULL;
    size_t n_a2freqs = 0;

    // ------------------------------------------------
    // IBDSEQ method implementation
    double ibd_errormax = 0.001;
    double ibd_errorprop = 0.25;
    double ibd_ibdlod = 3.0;
    double ibd_ibdtrim = 0.0;
    double ibd_max_error_array[5] = { 0.0 };
    //TODO implement
    uint64_t ibd_segment_max_n_missing_sites = 10000;

    double ibd_alpha = -1.0;
    double ibd_beta = -1.0;
    int ibd_dynamic_alpha = 0;
    int ibd_dynamic_beta = 0;

    // ------------------------------------------------
    // -> pair filters
    int allow_mispairs = -1;
    int pair_min_n_sites = -1;
    int min_n_pairs = -1;
    // ------------------------------------------------


    int seed = -1;
    int nThreads = 1;

    int rmInvarSites = 0;
    int rmMultiallelicSites = 0;

    int windowSize = 0;

    // input file type
    uint8_t in_ft = 0;

    kstring_t outfiles_list = KS_INITIALIZE;

};

//
/// @brief argStruct_get - parse command line arguments
/// @return pointer to the argStruct structure
argStruct* argStruct_get(int argc, char** argv);
void argStruct_destroy(argStruct* arg);


// print generic usage information
void print_help(FILE* fp);

/// print formula usage information; to be used in formula specific errors
void print_help_formula(FILE* fp);
#endif  // __ARG_STRUCT__
