#ifndef __VCF_READER__
#define __VCF_READER__

#include <htslib/kstring.h>
#include <htslib/tbx.h>

#include "argStruct.h"
#include "bootstrap.h"
#include "dataStructs.h"
#include "io.h"
#include "paramStruct.h"

/* FORWARD DECLARATIONS ----------------------------------------------------- */

typedef struct vcfData vcfData;
typedef struct gtData gtData;
typedef struct glData glData;

struct bootstrapReplicate;
struct blobStruct;
/* -------------------------------------------------------------------------- */

extern const int get_3x3_idx[3][3];

extern const int bcf_allele_charToInt[256];

int bcf_alleles_get_gtidx(int a1, int a2);
int bcf_alleles_get_gtidx(char a1, char a2);

struct vcfData {
    vcfFile *in_fp = NULL;
    bcf1_t *bcf = NULL;

    bcf_hdr_t *hdr = NULL;

    hts_idx_t *idx = NULL;
    hts_itr_t *itr = NULL;

    // tbx_t *tbx = NULL;

    int nContigs = 0;
    int nInd = 0;
    int nIndCmb = 0;
    int nSites = 0;
    int totSites = 0;

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
    int **JointGenoCountDistGT = NULL;

    // \def jgcd_gt[nBlocks][nIndCmb][nJointClasses]
    //      jgcd_gt[b][i][j] == number of sites where the ith pair of individuals have the jth joint genotype class in block b
    int ***jgcd_gt = NULL;

    // \def pair_shared_nSites[nBlocks][nIndCmb]
    //      pair_shared_nSites[b][i] == number of sites shared between the individuals in the ith pair in block b
    int **pair_shared_nSites = NULL;

    double **JointGenoCountDistGL = NULL;
    double **JointGenoProbDistGL = NULL;
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

    void init_JointGenoCountDistGL();
    void init_JointGenoProbDistGL();
    void init_JointGenoCountDistGT();

    void print_JointGenoCountDist();
    void print_JointGenoProbDist();

    void lngl_init(int doEM);
    void lngl_expand();

    void _print(FILE *fp);
    void _print();

    int records_next();
    int records_next(hts_itr_t *block_itr);

    void set_n_joint_categories(paramStruct *pars);
};

vcfData *vcfData_init(paramStruct *pars);
void vcfData_destroy(vcfData *v);

struct gtData {
    int32_t *data = NULL;

    int size_e = 0;
    int n_values = 0;
    int n_missing_ind = 0;

    // TODO nvalues vs ploidy
    int ploidy = 0;

    gtData(vcfData *vcfd);
    ~gtData();

    /// @brief ind pointer
    /// @param ind_i index of the individual
    /// @return  pointer to the start of the GTs for the given individual
    int *ind_ptr(const int ind_i);
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

    // TODO nvalues vs ploidy
    int ploidy = 0;
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

void readSites_GL(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt);
void readSites_GL(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt, blobStruct *blob);

void readSites_GT(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt);
void readSites_GT(vcfData *vcfd, paramStruct *pars, pairStruct **pairSt, blobStruct *blob);

int site_read_GL(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars, pairStruct **pairs);

int get_JointGenoDist_GT(const int contig_i, const int site_i, vcfData *vcf, paramStruct *pars);
int get_JointGenoDist_GT(const int contig_i, const int site_i, vcfData *vcfd, paramStruct *pars, const int block_i);

int GLtoGT_1_JointGenoDist(vcfData *vcf, paramStruct *pars);

int parse_VCF_GL(paramStruct *pars, vcfFile *in_fp, bcf_hdr_t *hdr, bcf1_t *bcf, blobStruct *blobSt);

#endif  // __VCF_READER__
