/**
 * @file    euclid.h
 * @brief   header file for euclid.cpp
 * @details contains functions for checking for Euclidean properties of distance matrices and Euclidean distance transformations
 */
#ifndef __EUCLID_H__
#define __EUCLID_H__

 /* INCLUDES ----------------------------------------------------------------- */
#include "shared.h"
/* END-OF-INCLUDES ---------------------------------------------------------- */

/* FORWARD-DECLARATIONS ----------------------------------------------------- */
/* END-OF-FORWARD-DECLARATIONS ---------------------------------------------- */

/* MACROS ------------------------------------------------------------------- */
///* END-OF-MACROS ------------------------------------------------------------ */


///* TYPEDEF-STRUCTS ---------------------------------------------------------- */
///* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

///* FUNCTION-DECLARATIONS ---------------------------------------------------- */

/// @brief cailliez - apply the Cailliez correction to a distance matrix
/// @param distmat   - pointer to the distance matrix (input)
/// @param n_elems   - number of elements in the distance matrix
/// @param corrected - pointer to the corrected distance matrix (output)
/// @param tole      - tolerance threshold
/// @param cor_zero  - bool indicating whether to correct values values close to zero (with tole as the threshold for being close to zero)
/// @return void
void cailliez(double* distmat, const size_t n_elems, double* corrected, const double tole, const bool cor_zero);

void test_cailliez(void);

/// @brief matrix_is_euclidean - check if a matrix is Euclidean
/// @param m    - pointer to the matrix
/// @param n    - number of elements in the matrix
/// @param tole - tolerance threshold
/// @return bool - true if the matrix is Euclidean, false otherwise
bool matrix_is_euclidean(double* m, int n, double tole);

/* END-OF-FUNCTION-DECLARATIONS --------------------------------------------- */

#endif  // __EUCLID_H__