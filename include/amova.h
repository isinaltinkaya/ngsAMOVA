/**
 * @file    amova.h
 * @brief   header file for amova.cpp
 * @details contains functions for performing Analysis of Molecular Variance (AMOVA)
 */
#ifndef __AMOVA_H__
#define __AMOVA_H__


 /* INCLUDES ----------------------------------------------------------------- */
#include "shared.h"
/* END-OF-INCLUDES ---------------------------------------------------------- */

/* FORWARD-DECLARATIONS ----------------------------------------------------- */
typedef struct metadata_t metadata_t;
typedef struct argStruct argStruct;
typedef struct dmat_t dmat_t;
typedef struct outfile_t outfile_t;

typedef struct amova_t amova_t;
typedef struct strArray strArray;

/* END-OF-FORWARD-DECLARATIONS ---------------------------------------------- */

/* MACROS ------------------------------------------------------------------- */
/* END-OF-MACROS ------------------------------------------------------------ */

/* TYPEDEF-STRUCTS ---------------------------------------------------------- */

/**
 * @brief amova - struct for storing AMOVA results
 * @note  if isShared==FALSE; then the arrays are allocated for each bootstrap replicate
 *        and the first array always stores the original data
 */
struct amova_t {

    /// ------------------------------------------///
    /// -> shared between bootstrap replicates <- ///

    /// @var metadata - pointer to metadata to be used in AMOVA
    metadata_t* metadata;

    /// @var  df  - array of degrees of freedom
    /// @note size = nLevels
    ///       df[i] = df for i-th level within the (i-1)-th level
    ///       [bootstrap] isShared = TRUE
    int* df;

    /// @var  df_total - total degrees of freedom
    ///       [bootstrap] isShared = TRUE
    int df_total;

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

    /// ------------------------------------------///
    /// -> unique to each bootstrap replicate <- ///

    size_t nRuns;

    /// @var  ss - array of sum of squares within (SS^(w))
    /// @note size = [nRuns][nLevels]
    ///       ss[i] = SS within the i-th level
    ///       ss[nLevels-1] = SS within individuals (currently unused! set to 0.0)
    ///       [bootstrap] isShared = FALSE
    double** ss;

    /// @var  ss_total - total sum of squares
    /// @note size = nRuns
    ///       [bootstrap] isShared = FALSE
    double* ss_total;

    /// @var  ssd - array of sum of squared deviations (SSD)
    /// @note size = [nRuns][nLevels]
    ///       [bootstrap] isShared = FALSE
    double** ssd;

    /// @var  ssd_total - total sum of squared deviations (SSD)
    /// @note size = nRuns
    ///       [bootstrap] isShared = FALSE
    double* ssd_total;

    /// @var  msd - array of mean squared deviations
    /// @note size = [nRuns][nLevels]
    ///       [bootstrap] isShared = FALSE
    double** msd;

    /// @var  msd_total - total mean squared deviation
    /// @note size = nRuns
    ///       [bootstrap] isShared = FALSE
    double* msd_total;

    /// @var  sigmasq - array of variance components
    /// @note size = [nRuns][nLevels]
    ///       sigmasq[r][i] = variance component for i-th level (for r-th run)
    ///       [bootstrap] isShared = FALSE
    double** sigmasq;

    /// @var  sigmasq_total - sum of all variance components
    /// @note size = nRuns
    ///       [bootstrap] isShared = FALSE
    double* sigmasq_total;

    /// @var  phi_xt - array of phi_xt statistics
    /// @note size = [nRuns][nLevels-1]
    ///       [bootstrap] isShared = FALSE
    double** phi_xt;

    /// @var  phi_xy - array of phi_xy statistics
    /// @note size = [nRuns][nLevels-2]
    ///       phi_xy[r][i] = Phi_{i(i-1)}^{*(r)} (r-th run)
    ///       [bootstrap] isShared = FALSE
    double** phi_xy;

    /// @note *_adj: allocated iff adjustment is needed (found a negative variance component)

    /// @var  sigmasq_adj - array of adjusted variance components 
    /// @note size = [nRuns][nLevels]
    ///       sigmasq_adj[r][i] = adjusted variance component for i-th level (for r-th run)
    ///       [bootstrap] isShared = FALSE
    double** sigmasq_adj;

    /// @var  sigmasq_total_adj - sum of all adjusted variance components
    /// @note size = [nRuns]
    double* sigmasq_total_adj;

    /// @var  phi_xt_adj - array of adjusted phi_xt statistics
    /// @note size = [nRuns][nLevels-1]
    ///       [bootstrap] isShared = FALSE
    double** phi_xt_adj;

    /// @var  phi_xy_adj - array of adjusted phi_xy statistics
    /// @note size = [nRuns][nLevels-2]
    ///       [bootstrap] isShared = FALSE
    double** phi_xy_adj;


};

/* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

/* FUNCTION-DECLARATIONS ----------------------------------------------------- */

amova_t* amova_get(dmat_t* dmat, metadata_t* mtd);
void amova_destroy(amova_t* amv);

/* END-OF-FUNCTION-DECLARATIONS ---------------------------------------------- */

// #endif  // __X_H__

/* ========================================================================== */


#endif  // __AMOVA_H__