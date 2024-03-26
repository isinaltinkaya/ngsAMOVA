/**
 * @file    amova.h
 * @brief   header file for amova.cpp
 * @details contains functions for performing Analysis of Molecular Variance (AMOVA)
 */
#ifndef __AMOVA_H__
#define __AMOVA_H__

#include <pthread.h>

#include "mathUtils.h"


typedef struct metadata_t metadata_t;
struct amovaBootstrapThreads;

typedef struct amova_t amova_t;
typedef struct amovaBootstrapThreads amovaBootstrapThreads;



/**
 * @brief amova - struct for storing AMOVA results
 * @note  if isShared==FALSE; then the arrays are allocated for each bootstrap replicate
 *        and the first array always stores the original data
 */
struct amova_t {

    metadata_t* metadata;

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

void amova_print_as_csv(amova_t* amv);


amova_t* amova_init(metadata_t* mtd, const int nAmovaRuns);

void amova_destroy(amova_t* amv);

amova_t* amova_get(dmat_t* dmat, metadata_t* mtd);



#endif  // __AMOVA_H__