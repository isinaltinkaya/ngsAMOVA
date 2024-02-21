#ifndef __AMOVA__
#define __AMOVA__

#include <pthread.h>

#include "mathUtils.h"


/// @brief GET_UPTRID_MATRIX_IJ - get the i,j element of an upper triangular matrix 
/// @note matrix is stored as a 1D array of size (n*(n+1))/2; including the diagonal
#define GET_UPTRID_MATRIX_IJ(i, j) ( ((i) > (j)) ? ( ((i)*((i)+1))/2 + (j) ) : ( ((j)*((j)+1))/2 + (i) ) )

struct distanceMatrixStruct;
typedef struct metadataStruct metadataStruct;
struct amovaBootstrapThreads;

typedef struct amovaStruct amovaStruct;
typedef struct amovaBootstrapThreads amovaBootstrapThreads;



/**
 * @brief amovaStruct - struct for storing AMOVA results
 * @note  if isShared==FALSE; then the arrays are allocated for each bootstrap replicate
 *        and the first array always stores the original data
 */
struct amovaStruct {

    metadataStruct* metadata;

    size_t nRuns;
    // size_t nPermutations; // TODO

    /// @var  df  - array of degrees of freedom
    /// @note size = nLevels
    ///       df[i] = df for i-th level within the (i-1)-th level
    ///       [bootstrap] isShared = TRUE
    int* df; // 

    /// @var  df_total - total degrees of freedom
    ///       [bootstrap] isShared = TRUE
    int df_total;

    /// @var  ss - array of sum of squares within (SS^(w))
    /// @note size = nLevels
    ///       ss[i] = SS within the i-th level
    ///       ss[nLevels-1] = SS within individuals (currently unused! set to 0.0)
    ///       [bootstrap] isShared = FALSE
    double** ss;

    /// @var  ss_total - total sum of squares
    ///       [bootstrap] isShared = FALSE
    double* ss_total;

    /// @var  ssd - array of sum of squared deviations (SSD)
    /// @note size = nLevels
    ///       [bootstrap] isShared = FALSE
    double** ssd;

    /// @var  ssd_total - total sum of squared deviations (SSD)
    ///       [bootstrap] isShared = FALSE
    double* ssd_total;

    /// @var  msd - array of mean squared deviations
    /// @note size = nLevels
    ///       [bootstrap] isShared = FALSE
    double** msd;

    /// @var  msd_total - total mean squared deviation
    ///       [bootstrap] isShared = FALSE
    double* msd_total;

    /// @var  vmat - used in the calculation of cmat
    /// @note size = (nLevels * (nLevels+1)) / 2
    ///       upper triangular matrix of v_ij elements used in calculation of variance coefficients
    ///       [bootstrap] isShared = TRUE
    double* vmat;

    /// @var  cmat - upper triangular matrix of variance coefficients
    /// @note size = (nLevels * (nLevels+1)) / 2
    ///       in classical amova permutation test, this is isShared = FALSE
    ///       in bootstrapping, this is isShared = TRUE since ninds are not permuted
    ///       [bootstrap] isShared = TRUE
    double* cmat;

    /// @var  lmat - inverse of cmat
    /// @note size = (nLevels * (nLevels+1)) / 2
    ///       [bootstrap] isShared = TRUE
    double* lmat;

    /// @var  sigmasq - array of variance components
    /// @note size = nLevels
    ///       sigmasq[i] = variance component for i-th level
    ///       [bootstrap] isShared = FALSE
    double** sigmasq;

    /// @var  sigmasq_total - sum of all variance components
    ///       [bootstrap] isShared = FALSE
    double* sigmasq_total;

    /// @var  phi_xt - array of phi_xt statistics
    /// @note size = nLevels-1
    ///       [bootstrap] isShared = FALSE
    double** phi_xt;

    /// @var  phi_xy - array of phi_xy statistics
    /// @note size = nLevels-2
    ///       phi_xy[i] = Phi_{i(i-1)}
    ///       [bootstrap] isShared = FALSE
    double** phi_xy;

};

void amovaStruct_print_as_csv(amovaStruct* amv);


inline amovaStruct* amovaStruct_init(metadataStruct* mtd, const int nAmovaRuns) {

    amovaStruct* ret = (amovaStruct*)malloc(sizeof(amovaStruct));
    ret->metadata = mtd;

    const size_t nLevels = (size_t)mtd->nLevels;

    const size_t nRuns = (size_t)nAmovaRuns;
    ret->nRuns = nRuns;


    const size_t nCmat = (nLevels * (nLevels + 1)) / 2;
    ret->vmat = NULL;
    ret->vmat = (double*)malloc((nCmat) * sizeof(double));
    ASSERT(ret->vmat != NULL);
    ret->cmat = NULL;
    ret->cmat = (double*)malloc((nCmat) * sizeof(double));
    ASSERT(ret->cmat != NULL);
    ret->lmat = NULL;
    ret->lmat = (double*)malloc((nCmat) * sizeof(double));
    ASSERT(ret->cmat != NULL);
    for (size_t i = 0; i < nCmat; ++i) {
        ret->cmat[i] = 0.0;
        ret->lmat[i] = 0.0;
        ret->vmat[i] = 0.0;
    }


    ret->df = NULL;
    ret->df = (int*)malloc((nLevels) * sizeof(int));
    ASSERT(ret->df != NULL);

    ret->df_total = 0;

    ret->ss = NULL;
    ret->ss = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->ss != NULL);

    ret->ss_total = NULL;
    ret->ss_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(ret->ss_total != NULL);

    ret->ssd = NULL;
    ret->ssd = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->ssd != NULL);

    ret->ssd_total = NULL;
    ret->ssd_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(ret->ssd_total != NULL);

    ret->msd = NULL;
    ret->msd = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->msd != NULL);

    ret->msd_total = NULL;
    ret->msd_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(ret->msd_total != NULL);

    ret->sigmasq = NULL;
    ret->sigmasq = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->sigmasq != NULL);

    ret->sigmasq_total = NULL;
    ret->sigmasq_total = (double*)malloc((nRuns) * sizeof(double));
    ASSERT(ret->sigmasq_total != NULL);

    ret->phi_xt = NULL;
    ret->phi_xt = (double**)malloc((nRuns) * sizeof(double*));
    ASSERT(ret->phi_xt != NULL);


    ret->phi_xy = NULL;
    if (nLevels > 2) {
        ret->phi_xy = (double**)malloc((nRuns) * sizeof(double*));
        ASSERT(ret->phi_xy != NULL);
    }

    for (size_t i = 0;i < nRuns;++i) {
        ret->ss[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(ret->ss[i] != NULL);

        ret->ss_total[i] = 0.0;

        ret->ssd[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(ret->ssd[i] != NULL);

        ret->ssd_total[i] = 0.0;

        ret->msd[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(ret->msd[i] != NULL);

        ret->msd_total[i] = 0.0;

        ret->sigmasq[i] = (double*)malloc((nLevels) * sizeof(double));
        ASSERT(ret->sigmasq[i] != NULL);

        ret->sigmasq_total[i] = 0.0;

        ret->phi_xt[i] = (double*)malloc((nLevels - 1) * sizeof(double));
        ASSERT(ret->phi_xt[i] != NULL);

        if (ret->phi_xy != NULL) {
            ret->phi_xy[i] = NULL;
            ret->phi_xy[i] = (double*)malloc((nLevels - 2) * sizeof(double));
            ASSERT(ret->phi_xy[i] != NULL);
            for (size_t j = 0;j < nLevels - 2;++j) {
                ret->phi_xy[i][j] = 0.0;
            }
        }


        for (size_t j = 0;j < nLevels;++j) {
            ret->ss[i][j] = 0.0;
            ret->ssd[i][j] = 0.0;
            ret->msd[i][j] = 0.0;
            ret->sigmasq[i][j] = 0.0;
            if (j < nLevels - 1) {
                ret->phi_xt[i][j] = 0.0;
            }
        }

    }

    return(ret);
}


inline void amovaStruct_destroy(amovaStruct* amv) {

    amv->metadata = NULL;

    for (size_t i = 0;i < amv->nRuns;++i) {
        FREE(amv->ss[i]);
        FREE(amv->ssd[i]);
        FREE(amv->msd[i]);
        FREE(amv->sigmasq[i]);
        FREE(amv->phi_xt[i]);
        if (amv->phi_xy != NULL) {
            FREE(amv->phi_xy[i]);
        }
    }

    FREE(amv->df);
    FREE(amv->ss);
    FREE(amv->ss_total);
    FREE(amv->ssd);
    FREE(amv->ssd_total);
    FREE(amv->msd);
    FREE(amv->msd_total);
    FREE(amv->vmat);
    FREE(amv->cmat);
    FREE(amv->lmat);
    FREE(amv->sigmasq);
    FREE(amv->sigmasq_total);
    FREE(amv->phi_xt);
    if (amv->phi_xy != NULL) {
        FREE(amv->phi_xy);
    }

    FREE(amv);

}


amovaStruct* amovaStruct_get(distanceMatrixStruct* dm, metadataStruct* mtd, blobStruct* blob);



#endif  // __AMOVA__