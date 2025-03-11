/**
 * @file    metadata.h
 * @brief   header file for metadata.cpp
 * @details contains metadata data structure (metadata_t) and associated functions
 */
#ifndef __METADATA_H__
#define __METADATA_H__


 /* INCLUDES ----------------------------------------------------------------- */
 /* END-OF-INCLUDES ---------------------------------------------------------- */

 /* FORWARD-DECLARATIONS ----------------------------------------------------- */
typedef struct metadata_t metadata_t;
typedef struct paramStruct paramStruct;
typedef struct strArray strArray;
typedef struct size_tArray size_tArray;

/* END-OF-FORWARD-DECLARATIONS ---------------------------------------------- */

/* MACROS ------------------------------------------------------------------- */
/* END-OF-MACROS ------------------------------------------------------------ */

/* TYPEDEF-STRUCTS ---------------------------------------------------------- */
/// @note
/// each new group name is added to groupNames and is given a group identifier index 
///
/// number of groups at level i = level2groupIndices[i]->len
/// index of the i-th group at level j = level2groupIndices[j]->d[i]
///
/// number of pairs of individuals that belong to group with index i = group2pairIndices[i]->len
/// array of pairs of individuals (pair indices) that belong to the group with index i = group2pairIndices[i]->d
///
/// number of individuals that belong to group with index i = group2indIndices[i]->len
/// array of individuals that belong to the group with index i = group2indIndices[i]->d
///
/// number of groups that are subgroups of the group with index i = group2subgroupIndices[i]->len
/// array of groups that are subgroups of the group with index i = group2subgroupIndices[i]->d
///
struct metadata_t {

    // total number of levels
    // e.g. individual, region, population, subpopulation
    //      nLevels = 4
    int nLevels;

    int nGroups;

    int nInd;

    // \def levelNames->d[nLevels]
    // names of levels in the hierarchy (e.g. region, population, subpopulation, individual)
    // in the hierarchical order
    // e.g. formula: "Individual ~ Region/Population/Subpopulation"
    //      levelNames = {"Region","Population","Subpopulation","Individual"}
    strArray* levelNames;

    // \def indNames->d[nInds]
    // indNames->d[i] = name of individual i (char *)
    /// @note free iff metadata==NULL && PROGRAM_HAS_INPUT_VCF ; else it is ptr to metadata->indNames
    strArray* indNames;

    // \def groupNames->d[nGroups]
    // N.B. does not contain individual names
    strArray* groupNames;

    size_tArray** level2groupIndices;

    size_tArray* group2levelIndices;

    size_tArray** group2indIndices;

    size_tArray** group2pairIndices;

    size_tArray** group2subgroupIndices;

};

/* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

/* FUNCTION-DECLARATIONS ----------------------------------------------------- */

metadata_t* metadata_init(const int in_nLevels);
metadata_t* metadata_read(paramStruct* pars);
void metadata_destroy(metadata_t* mtd);

/* END-OF-FUNCTION-DECLARATIONS ---------------------------------------------- */

#endif  // __METADATA_H__