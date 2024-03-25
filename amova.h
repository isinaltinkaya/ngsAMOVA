#ifndef __AMOVA__
#define __AMOVA__

#include <pthread.h>

#include "mathUtils.h"


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


amovaStruct* amovaStruct_init(metadataStruct* mtd, const int nAmovaRuns);

void amovaStruct_destroy(amovaStruct* amv);

amovaStruct* amovaStruct_get(paramStruct* pars, metadataStruct* mtd);



#endif  // __AMOVA__