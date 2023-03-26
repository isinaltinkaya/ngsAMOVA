#ifndef __DXY__
#define __DXY__

#include "mathUtils.h"
#include "dataStructs.h"
#include "io.h"









/// @brief estimate_dxy_2groups - estimate pairwise dxy between two groups
// double estimate_dxy_2groups(const int glob_idx1, const int glob_idx2, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars);


// Dxy analysis entry point  
void doDxy(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd);


typedef struct dxyStruct
{
    kstring_t *kbuf;

    // number of dxy values
    int nDxy=0;
    
    size_t _dxy=100; // initial value for malloc; will be increased if needed

    double* dxy;
    char** groupNames1;
    char** groupNames2;
    char** levelNames;

	void print(IO::outputStruct *out_dxy_fs);

    void print_struct();

    dxyStruct(const int printDxy);
    ~dxyStruct();

    /// @brief estimate_dxy_2groups - estimate pairwise dxy between two groups
    /// @details edits kbuf inplace
    void estimate_dxy_2groups(const int local_idx1, const int local_idx2, const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars);

    /// @brief estimate_dxy_allGroupsAtLevel - estimate pairwise dxy between all groups at a given level
    void estimate_dxy_allGroupsAtLevel(const int lvl, distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars);


    /// @brief estimate_dxy_allLevels - estimate pairwise dxy between all groups at all levels
    void estimate_dxy_allLevels(distanceMatrixStruct *dMS, metadataStruct *mtd, paramStruct *pars);

    void expand();


} dxyStruct;

// read dxyStruct from dxy file
dxyStruct *dxyStruct_read(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd);
dxyStruct *dxyStruct_get(argStruct *args, paramStruct *pars, distanceMatrixStruct *dMS, metadataStruct *mtd);

#endif