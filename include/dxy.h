/**
 * @file    dxy.h
 * @brief   header file for dxy.cpp
 * @details contains dxy analysis related data structures and functions
 */
#ifndef __DXY_H__
#define __DXY_H__


/* INCLUDES ----------------------------------------------------------------- */
#include "shared.h"
/* END-OF-INCLUDES ---------------------------------------------------------- */

/* FORWARD-DECLARATIONS ----------------------------------------------------- */
typedef struct metadata_t metadata_t;
typedef struct dmat_t dmat_t;
/* END-OF-FORWARD-DECLARATIONS ---------------------------------------------- */

/* MACROS ------------------------------------------------------------------- */
/* END-OF-MACROS ------------------------------------------------------------ */

/* TYPEDEF-STRUCTS ---------------------------------------------------------- */
/* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

/* FUNCTION-DECLARATIONS ----------------------------------------------------- */
/* END-OF-FUNCTION-DECLARATIONS ---------------------------------------------- */


//TODO use dmats only, rm dxy data d labels etc
// ur not using the dmats right now 
//TODO add test 

/// @brief dxy_t - structure for dxy analysis results
///
/// @field kbuf - kstring buffer for printing
/// @field n - number of dxy values
/// @field dxy - array of dxy values
/// @field g1names_p - array of group names for group 1 corresponding to dxy values
/// @field g2names_p - array of group names for group 2 corresponding to dxy values
/// @field levelnames_p - array of level names corresponding to dxy values
///
/// @details
/// dxy[i] is the dxy value between g1names_p[i] and g2names_p[i] (at levelnames_p[i])
///
typedef struct dxy_t {

    /// @var nLevels - number of levels in the metadata (excluding the "individual" level)
    size_t nLevels;

    /// @var n - number of dxy values (number of pairwise group comparisons)
    size_t n;


    /// @var d - dxy data (array of dxy values)
    double* d;

    /// @var dmat - array of distance matrices per level
    /// @size mtd->nLevels - 1 (excluding the "individual" level)
    /// @note dmat[i] can be NULL if level is skipped (e.g. only one group at that level)
    /// @note check if dmat[i] is NULL before using/freeing it
    dmat_t** dmat;

    /// @var g1names_p - array of pointers to group names for group 1 corresponding to the dxy values
    /// @size n
    /// @note only the pointers are stored, the actual group names are stored in the metadata_t
    /// @note free only the array, not the pointers
    char** g1names_p;

    /// @var g2names_p - array of pointers to group names for group 2 corresponding to the dxy values
    /// @size n
    /// @note only the pointers are stored, the actual group names are stored in the metadata_t
    /// @note free only the array, not the pointers
    char** g2names_p;


    /// @var levelnames_p - array of pointers to level names corresponding to the dxy values
    /// @size n
    /// @note only the pointers are stored, the actual level names are stored in the metadata_t
    /// @note free only the array, not the pointers
    char** levelnames_p;

} dxy_t;

dxy_t* dxy_read(paramStruct* pars, dmat_t* dmat, metadata_t* mtd);
dxy_t* dxy_get(dmat_t* dmat, metadata_t* mtd);

void dxy_destroy(dxy_t* dxy);




#endif // __DXY_H__