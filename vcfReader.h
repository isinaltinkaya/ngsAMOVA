#ifndef __VCF_READER__
#define __VCF_READER__

#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <stdbool.h>

#include "bootstrap.h"
#include "io.h"
#include "paramStruct.h"

// UNMOD: unmodified, do not perform any transform or normalization
#define ARG_GLDATA_TYPE_UNMOD 0
#define ARG_GLDATA_TYPE_GL 1
#define ARG_GLDATA_TYPE_LNGL 2
#define ARG_GLDATA_TYPE_LOG10GL 3


// BASES are: { A, C, G, T, BASE_UNOBSERVED }
// base indices: A:0 C:1 G:2 T:3 BASE_UNOBSERVED:4
// where BASE_UNOBSERVED is the unobserved allele denoted by <*> or <NON_REF>
#define BASE_UNOBSERVED 4


/* FORWARD DECLARATIONS ----------------------------------------------------- */

typedef struct vcfData vcfData;
typedef struct gtData gtData;

struct blobStruct;
typedef struct bblocks_t bblocks_t;
/* -------------------------------------------------------------------------- */


/// @brief nDerToM33Idx - convert the number of derived alleles to the index of the corresponding combination in the 3x3 matrix
/// @details
//
// 0 1 2
// 00 01 02
// MMMM MMMm MMmm
//
// 3 4 5
// 10 11 12
// MmMM MmMm Mmmm
//
// 6 7 8
// 20 21 22
// mmMM mmMm mmmm
extern const int nDerToM33Idx[3][3];

/// @brief acgt_charToInt - convert a char base to int acgt index (internal representation)
/// @details
/// bases are internally represented as integers
///
/// a,c,g,t,n
/// A,C,G,T,N
/// 0,1,2,3,4
extern const int acgt_charToInt[256];

int bcf_alleles_get_gtidx(int a1, int a2);
int bcf_alleles_get_gtidx(char a1, char a2);

// index file types
#define IDX_NONE 0
#define IDX_CSI (1 << 0)  // 1
#define IDX_TBI (1 << 1)  // 2

/// @brief require_index - check if the analysis requires loading an index file
/// @param pars
/// @return 0 if no index file is required, otherwise return an IDX_* value
///         indicating which index file type is required given the input file type
int require_index(paramStruct* pars);

/// @brief require_unpack - check if the analysis requires unpacking
/// @return 0 if no unpacking is required, otherwise return a BCF_UN_* value
int require_unpack();


typedef struct gldata_t gldata_t;
struct gldata_t {

    /// @details known a priori (args->doEM)
    uint8_t n_gls;

    /// @details known a priori (pars->nInd; from vcfd->hdr->n[BCF_DT_SAMPLE])
    size_t n_ind;

    /// @details at runtime, worst case scenario: total number of sites in the vcf file header for all contigs+step_size
    /// at the end of read_sites, equal to pars->nSites
    size_t size;

    /// @var step_size - step size for reallocating/expanding the number of gl arrays
    /// step_size*n_gls is the added size when reallocating
    size_t step_size;

    /// @var type - type of genotype likelihoods
    /// @details values: ARG_GLDATA_TYPE_*
    uint8_t type;

    /// @var d[n_ind][n_sites][n_gls]
    /// @details
    /// d[i] - (double**) array of all gls for individual i 
    /// d[i][j] - (double*) array of gls for individual i at site j
    ///           allocated iff !mis[i][j]; size=n_gls
    ///           otherwise NULL
    /// d[i][j][k] - (double) gl for individual i at site j for genotype k
    double*** d;

    /// @var mis[n_ind][n_sites] - indicator for whether if an individual has missing data at a site
    /// @details
    /// mis[i] - array of bool indicators for missing data for individual i
    /// mis[i][j] - bool indicator for missing data for individual i at site j
    /// @note 
    /// mis[i][j] == true if individual i does not have data at site j, 
    /// OR, if it did have data but due to some filter we decided to ignore it
    bool** mis;

};

struct vcfData {

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

    gldata_t* gldata = NULL;

    int nGT = 0;


    // index of the unobserved allele in vcfd->rec->d.allele (if any)
    // set to -1 if no unobserved allele is found in d.alleles
    int allele_unobserved = -1;

    int nJointClasses = 0;

    void _print(FILE* fp);
    void _print();

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

vcfData* vcfData_init(paramStruct* pars, metadataStruct* metadata);
void vcfData_destroy(vcfData* v);


/*
 * @template struct get_data
 *
 * @abstract		wrapper for bcf_get_data_*
 *
 * @field data
 *
 * @field size_e	watermark for number of elements
 *
 * @field n			number of returned values
 * 					if <0; error
 * 					else; number of written values
 * 					used for accessing entries
 *
 */
template <typename T>
struct get_data {
    T* data = NULL;

    int size_e = 0;
    int n = 0;
    int n_missing_ind = 0;

    int ploidy = 0;  // ploidy = n_values / nInd
    int n_values = 0;

    T& operator[](unsigned i) {
        return data[i];
    }

    bool is_empty() const {
        return data == NULL;
    }

    ~get_data() {
        FREE(data);
    }
};

void readSites(jgtmat_t* jgtm, bblocks_t* bblocks, vcfData* vcfd, paramStruct* pars);



int read_site_with_alleles(const int site_i, vcfData* vcfd, paramStruct* pars);

#endif  // __VCF_READER__
