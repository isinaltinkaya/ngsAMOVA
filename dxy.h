#ifndef __DXY__
#define __DXY__

#include "dataStructs.h"
#include "io.h"
#include "mathUtils.h"



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

    /// @var dm - array of distance matrices per level
    /// @size mtd->nLevels - 1 (excluding the "individual" level)
    /// @note dm[i] can be NULL if level is skipped (e.g. only one group at that level)
    /// @note check if dm[i] is NULL before using/freeing it
    dmat_t** dm;

    /// @var g1names_p - array of pointers to group names for group 1 corresponding to the dxy values
    /// @size n
    /// @note only the pointers are stored, the actual group names are stored in the metadataStruct
    /// @note free only the array, not the pointers
    char** g1names_p;

    /// @var g2names_p - array of pointers to group names for group 2 corresponding to the dxy values
    /// @size n
    /// @note only the pointers are stored, the actual group names are stored in the metadataStruct
    /// @note free only the array, not the pointers
    char** g2names_p;


    /// @var levelnames_p - array of pointers to level names corresponding to the dxy values
    /// @size n
    /// @note only the pointers are stored, the actual level names are stored in the metadataStruct
    /// @note free only the array, not the pointers
    char** levelnames_p;

} dxy_t;

dxy_t* dxy_read(paramStruct* pars, dmat_t* dmat, metadataStruct* mtd);
dxy_t* dxy_get(paramStruct* pars, dmat_t* dmat, metadataStruct* mtd);

inline void dxy_destroy(dxy_t* dxy) {
    FREE(dxy->d);
    FREE(dxy->g1names_p);
    FREE(dxy->g2names_p);
    FREE(dxy->levelnames_p);
    for (size_t i = 0;i < dxy->nLevels;++i) {
        if (dxy->dm[i] != NULL) {
            dmat_destroy(dxy->dm[i]);
        }
    }
    FREE(dxy->dm);
    FREE(dxy);
    return;
}

void dxy_print(dxy_t* dxy);



#endif // __DXY__