#ifndef __VCF_READER__
#define __VCF_READER__

#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <stdbool.h>

#include "bootstrap.h"
#include "io.h"
#include "paramStruct.h"

/* FORWARD DECLARATIONS ----------------------------------------------------- */

typedef struct vcfData vcfData;
typedef struct gtData gtData;
typedef struct glData glData;

struct blobStruct;
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
int require_index(paramStruct *pars);

/// @brief require_unpack - check if the analysis requires unpacking
/// @return 0 if no unpacking is required, otherwise return a BCF_UN_* value
int require_unpack();

void set_contig_vcfToAlleles(vcfData *vcfd, int **contig_vcfToAlleles, allelesStruct *alleles);

struct vcfData {
    vcfFile *in_fp = NULL;
    bcf1_t *rec = NULL;

    bcf_hdr_t *hdr = NULL;

    hts_idx_t *idx = NULL;
    hts_itr_t *itr = NULL;

    int require_index = 0;

    int nContigs = 0;
    int *skipContigs = NULL;

    int nInd = 0;
    int nIndCmb = 0;

    int *contig_vcfToAncDer = NULL;
    int *contig_vcfToMajMin = NULL;

    /// @brief get_contigIdx_in_alleles - get the index of a given contig in the allelesStruct
    /// @param contig_vcfToAlleles
    /// @param alleles
    /// @return index of the contig_i in the allelesStruct if found, otherwise -1
    const int get_contigIdx_in_alleles(int **contig_vcfToAlleles, allelesStruct *alleles);

    /*
     * lngl[nSites][nInd*nGT*double]
     * genotype likelihoods in natural log
     * for each individual at each site
     *
     * nGT=10 == store all 10 values per individual
     * nGT=3 == store only 3 values per individual
     * 		corresponding to (0,0), (0,1), (1,1)
     */
    double **lngl = NULL;
    size_t _lngl = 4096;
    int nGT = 0;

    /*
     * [nIndCmb][9+1]
     * last element contains total number of sites shared
     */
    int **JointGenotypeCountMatrixGT = NULL;

    // \def jgcd_gt[nBlocks][nIndCmb][nJointClasses]
    //      jgcd_gt[b][i][j] == number of sites where the ith pair of individuals have the jth joint genotype class in block b
    int ***jgcd_gt = NULL;

    // \def pair_shared_nSites[nBlocks][nIndCmb]
    //      pair_shared_nSites[b][i] == number of sites shared between the individuals in the ith pair in block b
    int **pair_shared_nSites = NULL;

    double **JointGenotypeCountMatrixGL = NULL;
    double **JointGenotypeProbMatrixGL = NULL;
    int nJointClasses = 0;

    // TODO instead of copying the names, just store the order
    // and access the names via bcf_hdr_id2name(hdr, i) where i is the order of the individual
    // in the bcf file

    char **indNames = NULL;  // individual names in bcf order

    void addIndNames() {
        indNames = (char **)malloc(nInd * sizeof(char *));
        for (size_t i = 0; i < (size_t)nInd; i++) {
            indNames[i] = NULL;
            indNames[i] = strdup(hdr->samples[i]);
        }
    }

    // given number of genotype categories (nGT_)
    // set nGT to given nGT, set nJointClasses to nGT*nGT
    // void set_nGT(const int nGT_);

    void init_JointGenotypeCountMatrixGL();
    void init_JointGenotypeProbMatrixGL();
    void init_JointGenotypeCountMatrixGT();

    void print_JointGenotypeCountMatrix();

    void lngl_init(int doEM);
    void lngl_expand();

    void _print(FILE *fp);
    void _print();

    /// @brief records_next - get the next record
    /// @return 0 if there are no more records, otherwise return 1
    int records_next();

    void set_n_joint_categories(paramStruct *pars);

    /// @brief unpack - unpack the bcf record based on the value from require_unpack()
    void unpack(void);

    gtData *gts = NULL;
    glData *gls = NULL;

    void site_gts_get(const int a1, const int a2);
    void site_get_gls(void);

    /// @brief get_rec_contig_id - get the contig id of the current record
    const char *get_contig_name(void);

    /// @brief get_rec_contig_id(i) - get the contig id of the contig with id i
    const char *get_contig_name(const int32_t i);
};

vcfData *vcfData_init(paramStruct *pars);
void vcfData_destroy(vcfData *v);

struct gtData {
    int32_t *data = NULL;

    int size_e = 0;
    int n_values = 0;
    int n_missing_ind = 0;

    int ploidy = 0;  // ploidy = n_values / nInd

    gtData();
    ~gtData();

    /// @brief ind pointer
    /// @param ind_i index of the individual
    /// @return  pointer to the start of the GTs for the given individual
    int *ind_ptr(const int ind_i);

    /// @brief pass_minInd_threshold - check if the minInd threshold is met
    /// @param nInd number of individuals
    /// @return true if the minInd threshold is met, false otherwise
    bool pass_minInd_threshold(const int nInd);

    /// @brief get_n_derived_alleles_ind - get the number of derived alleles for the given individual
    /// @param ind_i index of the individual
    /// @return number of derived alleles found in the genotype of the given individual
    /// -1 if the individual is missing
    int get_n_derived_alleles_ind(const int ind_i, char **site_alleles);
    int get_n_derived_alleles_ind(const int ind_i);

    // internal representation: ACGT -> 0123
    // bcf representation : RefAllele,AltAllele1,AltAllele2,AltAllele3 -> xxxx

    // a1 == ancestral/major allele = 0
    // a2 == derived/minor allele = 1
    // \def intbase2state[4] == lookup table for converting the internal acgt order to allelic state, where 0 is ancestral, 1 is derived, -1 is other
    //      in the order {A, C, G, T}
    //      and initialized as {-1, -1, -1. -1}
    //      then alleles set to their corresponding allele type
    // @example if major is G and minor is C, then intbase2state = {1, 0, -1, -1}
    int *intbase2state = NULL;

    // \def acgt2alleles[4] == lookup table for converting the internal acgt index to the allele index in rec->d.allele
    //     in the order {A, C, G, T}
    //     and initialized as {-1, -1, -1, -1}
    //     then alleles set to their corresponding index in rec->d.allele
    // @example if at a site REF=A ALT=G,T, then acgt2alleles = {0, -1, 1, 2}

    // int *acgt2alleles = NULL;

    // \def alleles2intbase[4] == lookup table for converting the allele index in rec->d.allele to the internal acgt index

    // alleles is the allele index in bcf->d.allele(n=bcf->n_allele)

    // \def alleleidx2state[nAlleles] == lookup table for converting the allele index in rec->d.allele to the allelic state
    //      in the order of rec->d.allele
    //      and initialized as {-1, -1, -1, -1}
    //      then alleles set to their corresponding allele type
    // @example if major is G and minor is C, then alleleidx2state = {1, 0, -1, -1}
    int *alleleidx2state = NULL;

    /// @brief get_alleleidx2state - get the allelic state for the given allele index in rec->d->allele (alleleidx)
    /// @param alleleidx index of the allele in rec->d->allele
    /// @return the allelic state for the given allele index (0 for ancestral/major, 1 for derived/minor, otherwise -1)
    int get_alleleidx2state(const int alleleidx);

    void site_gts_clear();
};

/// @brief  genotype likelihood data
///
struct glData {
    float *data = NULL;

    int size_e = 0;
    int n_values = 0;
    int n_gls = 0;  // number of genotype likelihoods per individual (3 or 10)j
    int n_missing_ind = 0;

    glData(vcfData *vcfd);
    ~glData();

    int ind_data_isMissing(const int ind_i);

    /// @brief ind pointer
    /// @param ind_i index of the individual
    /// @return  pointer to the start of the GLs for the given individual
    float *ind_ptr(const int ind_i);
};

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
    T *data = NULL;

    int size_e = 0;
    int n = 0;
    int n_missing_ind = 0;

    int ploidy = 0;  // ploidy = n_values / nInd
    int n_values = 0;

    T &operator[](unsigned i) {
        return data[i];
    }

    bool is_empty() const {
        return data == NULL;
    }

    ~get_data() {
        FREE(data);
    }
};

void readSites(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt, blobStruct *blob);
void readSites(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt);

// @return
// 1    skip site
int site_read_GL(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars, pairStruct **pairs);

// return 1 if skipped
// block_i  -1 if block bootstrapping is disabled
int get_JointGenotypeMatrix_GT(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars, const int block_i);

int GLtoGT_1_JointGenotypeMatrix(vcfData *vcf, paramStruct *pars);

/// @brief read_site_with_alleles - for a given site fetch the alleles at site from allelesStruct
/// @return
/// 0    if positions match and everything works fine
/// 1    if 'contig:site' in VCF cannot be found in alleles file
/// 2    if a1 or a2 is not a single character
int read_site_with_alleles(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars, int *a1, int *a2);

#endif  // __VCF_READER__
