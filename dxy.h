#ifndef __DXY__
#define __DXY__

#include "dataStructs.h"
#include "io.h"
#include "mathUtils.h"

// TODO dxystruct is also a distancematrixstruct template

/// @brief dxyStruct - structure for dxy analysis results
///
/// @field kbuf - kstring buffer for printing
/// @field nDxy - number of dxy values
/// @field dxy - array of dxy values
/// @field groupNames1 - array of group names for group 1 corresponding to dxy values
/// @field groupNames2 - array of group names for group 2 corresponding to dxy values
/// @field levelNames - array of level names corresponding to dxy values
///
/// @details
/// dxy[i] is the dxy value between groupNames1[i] and groupNames2[i] (at levelNames[i])
///
typedef struct dxyStruct {
    // number of dxy values
    int nDxy = 0;

    size_t _dxyArr = 100;  // initial value for malloc; will be increased if needed

    double *dxyArr;
    char **groupNames1;
    char **groupNames2;
    char **levelNames;

    void print(IO::outputStruct *out_dxy_fs);
    void print_struct();
    void expand();

    dxyStruct();
    ~dxyStruct();

    /// @brief estimate_dxy_2groups - estimate pairwise dxy between two groups
    /// @return number of dxy values estimated (==1 since only one pair of groups:
    int estimate_dxy_2groups(const int local_idx1, const int local_idx2, const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars);

    /// @brief estimate_dxy_allGroupsAtLevel - estimate pairwise dxy between all groups at a given level
    /// @return number of dxy values estimated (==number of pairs of groups)
    int estimate_dxy_allGroupsAtLevel(const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars);

    /// @brief estimate_dxy_allLevels - estimate pairwise dxy between all groups at all levels
    /// @return number of dxy values estimated (==number of pairs of groups)
    int estimate_dxy_allLevels(distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars);

} dxyStruct;

// read dxyStruct from dxy file
dxyStruct *dxyStruct_read(paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd);
dxyStruct *dxyStruct_get(paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd);

#endif