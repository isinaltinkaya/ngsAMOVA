#ifndef __DXY__
#define __DXY__

#include "mathUtils.h"
#include "dataStructs.h"
#include "io.h"









/// @brief estimate_dxy_2groups - estimate pairwise dxy between two groups
// double estimate_dxy_2groups(const int glob_idx1, const int glob_idx2, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars);

/// @brief estimate_dxy_2groups - estimate pairwise dxy between two groups
/// @details edits kbuf inplace
void estimate_dxy_2groups(const int local_idx1, const int local_idx2, const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars, kstring_t *kbuf);

/// @brief estimate_dxy_allGroupsAtLevel - estimate pairwise dxy between all groups at a given level
void estimate_dxy_allGroupsAtLevel(const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars, kstring_t *kbuf);


/// @brief estimate_dxy_allLevels - estimate pairwise dxy between all groups at all levels
void estimate_dxy_allLevels(distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars, kstring_t *kbuf);


// Dxy analysis entry point  
void doDxy(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd);

#endif