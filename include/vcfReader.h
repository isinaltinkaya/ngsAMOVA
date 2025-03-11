/**
 * @file    vcfReader.h
 * @brief   header file for vcfReader.cpp
 * @details contains vcfReader related data structures and functions
 */

#ifndef __VCF_READER__
#define __VCF_READER__

/* INCLUDES ----------------------------------------------------------------- */

#include "io.h"
#include "shared.h"
#include "paramStruct.h"

/* END-OF-INCLUDES ---------------------------------------------------------- */


/// @brief acgt_charToInt - convert a char base to int acgt index (internal representation)
/// @details
/// bases are internally represented as integers
///
/// a,c,g,t,n
/// A,C,G,T,N
/// 0,1,2,3,4
extern const int acgt_charToInt[256];



/* FORWARD-DECLARATIONS ----------------------------------------------------- */

typedef struct vcfData vcfData;
typedef struct bblocks_t bblocks_t;
typedef struct ibds_t ibds_t;

/* END-OF-FORWARD-DECLARATIONS ------------------------------------------------ */




/* MACROS ------------------------------------------------------------------- */


// BASES are: { A, C, G, T, BASE_UNOBSERVED }
// base indices: A:0 C:1 G:2 T:3 BASE_UNOBSERVED:4
// where BASE_UNOBSERVED is the unobserved allele denoted by <*> or <NON_REF>
#define BASE_UNOBSERVED 4

// index file types
#define IDX_NONE    (0 << 0)  // 0
#define IDX_CSI     (1 << 0)  // 1
#define IDX_TBI     (1 << 1)  // 2

#define SITE_WILL_BE_SKIPPED(reason) (((reason) != SKIP_SITE_REASON_NOSKIP))

#define GLDATA_GET_GLSPTR_IND_SITE(gldata, ind, site) ((gldata)->d[(ind)] + ((site) * (gldata)->n_gls))

#define GLDATA_GET_GLPTR_IND_SITE_GL(gldata, ind, site, gl) ((gldata)->d[(ind)] + ((site) * (gldata)->n_gls + (gl)))

#define GLDATA_GET_INDPTR_AT_SITE(gldata, ind, site) ((gldata)->d[(ind)] + ((site) * (gldata)->n_gls))


/* END-OF-MACROS ------------------------------------------------------------- */

/* TYPEDEF-STRUCTS ---------------------------------------------------------- */

typedef struct gldata_t gldata_t;
struct gldata_t {

    /// @var n_ind - number of individuals
    size_t n_ind;

    /// @var n_sites - number of sites
    size_t n_sites;

    /// @var n_gls - number of genotype likelihoods
    size_t n_gls;

    /// @var step_size - step size in terms of n_sites for reallocating
    /// @details used for reallocating d[i] for i in [0, n_ind)
    size_t step_size;


    /// @var d - 2D array of genotype likelihoods
    /// @details
    // [memory layout]
    // data[n_ind] -> is a ptr to -> data[n_sites*n_gls]
    // for each individual i in [0, n_ind);
    // the data_i array is a linearized 2D array of [n_sites][n_gls]
    // storing the genotype likelihoods for each site as follows:
    // [{gl1s1}, ..., {glAs1}][{gl1s2}, ..., {glAs2}]...[{gl1sN}, ..., {glAsN}]
    // where A in glA is the number of genotype likelihoods
    // and N is the number of sites
    // and s is the site index
    float** d;

};

/* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

/* FUNCTIONS ---------------------------------------------------------------- */
static inline size_t gldata_estimate_max_mem_use(const size_t gldata_n_ind, const size_t gldata_n_sites, const size_t gldata_n_gls) {
    size_t size = sizeof(gldata_t);
    size += GET_ARR_SIZE_BYTES_2D(float, gldata_n_ind, gldata_n_sites * gldata_n_gls);
    return(size);
}


static inline gldata_t* gldata_init(const size_t n_ind, const size_t n_sites, const size_t n_gls, const size_t step_size) {

    // -> init
    gldata_t* gldata = (gldata_t*)malloc(sizeof(gldata_t));
    ASSERT(gldata != NULL);
    gldata->n_gls = 0;
    gldata->n_ind = 0;
    gldata->n_sites = 0;
    gldata->step_size = 0;
    gldata->d = NULL;

    // -> alloc/set
    gldata->n_gls = n_gls;
    gldata->n_ind = n_ind;
    gldata->n_sites = n_sites;
    gldata->step_size = step_size;

    const size_t size_each = n_sites * n_gls;
    gldata->d = (float**)malloc(n_ind * sizeof(float*));
    ASSERT(gldata->d != NULL);
    for (size_t i = 0; i < n_ind; ++i) {
        gldata->d[i] = (float*)malloc(size_each * sizeof(float));
        ASSERT(gldata->d[i] != NULL);
        for (size_t j = 0; j < size_each; ++j) {
            gldata->d[i][j] = 0.0; // initval
        }
    }

    return(gldata);

}

static inline void gldata_destroy(gldata_t* gldata) {
    for (size_t i = 0; i < gldata->n_ind; ++i) {
        FREE(gldata->d[i]);
    }
    FREE(gldata->d);
    FREE(gldata);
    return;
}


/// @brief gldata_expand - expand the per-ind data arrays n_sites by step_size
static inline void gldata_expand(gldata_t* gldata) {

    const size_t old_size_perind = gldata->n_sites * gldata->n_gls;
    const size_t new_size_perind = old_size_perind + (gldata->step_size * gldata->n_gls);

    for (size_t i = 0; i < gldata->n_ind; ++i) {
        gldata->d[i] = (float*)realloc(gldata->d[i], new_size_perind * sizeof(float));
        ASSERT(gldata->d[i] != NULL);
        for (size_t j = old_size_perind; j < new_size_perind; ++j) {
            gldata->d[i][j] = 0.0; // initval
        }
    }
    gldata->n_sites += gldata->step_size;
    return;
}

struct vcfData {

    size_t max_nsites = 0;
    size_t max_percontig_nsites = 0;

    vcfFile* in_fp = NULL;

    bcf1_t* rec = NULL;
    bcf_hdr_t* hdr = NULL;

    hts_idx_t* idx = NULL;
    tbx_t* tbx = NULL;

    hts_itr_t* itr = NULL;

    // @param DO_BCF_UNPACK
    // determines the level of unpacking needed for the specified analysis
    // init at vcfData_init()
    // used at every rec read
    // @values one of BCF_UN_* values
    // BCF_UN_STR: up to ALT inclusive
    // BCF_UN_FLT: up to FILTER
    // BCF_UN_INFO: up to INFO
    // BCF_UN_SHR: all shared information
    // BCF_UN_FMT: unpack format and each sample
    // BCF_UN_ALL: everything
    int DO_BCF_UNPACK;

    int nContigs = 0;
    int* skipContigs = NULL;


    // index of the unobserved allele in vcfd->rec->d.allele (if any)
    // set to -1 if no unobserved allele is found in d.alleles
    int allele_unobserved = -1;

    /// @brief records_next - get the next record
    /// @return 0 if there are no more records, otherwise return 1
    int records_next();


    /// @brief unpack - unpack the bcf record based on the value from require_unpack()
    void unpack(void);

    /// @brief get_rec_contig_id - get the contig id of the current record
    const char* get_contig_name(void);

    /// @brief get_rec_contig_id(i) - get the contig id of the contig with id i
    const char* get_contig_name(const int32_t i);

};

vcfData* vcfData_init(paramStruct* pars, metadata_t* metadata);
void vcfData_destroy(vcfData* v);

void readSites(vcfData* vcfd, paramStruct* pars, jgtmat_t* jgtmat, bblocks_t* bblocks, ibds_t* ibds, gldata_t* gldata);

/// @brief require_index - check if the analysis requires loading an index file
/// @param pars
/// @return 0 if no index file is required, otherwise return an IDX_* value
///         indicating which index file type is required given the input file type
int require_index(paramStruct* pars);

/// @brief require_unpack - check if the analysis requires unpacking
/// @return 0 if no unpacking is required, otherwise return a BCF_UN_* value
int require_unpack();

#endif  // __VCF_READER__