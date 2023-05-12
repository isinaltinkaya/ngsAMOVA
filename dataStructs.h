#ifndef __DATA_STRUCTS__
#define __DATA_STRUCTS__

#include "io.h"
#include "mathUtils.h"
#include "paramStruct.h"
#include "vcfReader.h"

/* FORWARD DECLARATIONS ----------------------------------------------------- */

typedef struct pairStruct pairStruct;
typedef struct distanceMatrixStruct distanceMatrixStruct;
typedef struct metadataStruct metadataStruct;
typedef struct vcfData vcfData;
typedef struct indPairThreads indPairThreads;

void spawnThreads_pairEM(argStruct *args, paramStruct *pars, pairStruct **pairSt, vcfData *vcfd, distanceMatrixStruct *distMatrix);
void setInputFileType(paramStruct *pars, int inFileType);
/* -------------------------------------------------------------------------- */

// prepare distance matrix using genotype likelihoods
void get_distanceMatrix_GL(argStruct *args, paramStruct *pars, distanceMatrixStruct *distanceMatrix, vcfData *vcfd, pairStruct **pairSt);
// prepare distance matrix using genotypes
void get_distanceMatrix_GT(argStruct *args, paramStruct *pars, distanceMatrixStruct *distanceMatrix, vcfData *vcfd, pairStruct **pairSt);

/// trim spaces from the beginning and end of a char* (inplace)
/// @param str - char* to trim
void trimSpaces(char *str);

struct pairStruct {
    int idx;  // index of the pair
    int i1;   // index of the first individual
    int i2;   // index of the second individual

    // number of shared sites
    size_t snSites = 0;

    // contains index of the shared sites of the pair
    int *sharedSites = NULL;
    size_t _sharedSites = 4096;

    double d;
    int n_em_iter;

    double *optim_jointGenoCountDist = NULL;
    double *optim_jointGenoProbDist = NULL;

    pairStruct(paramStruct *pars_, int pidx, int i1_, int i2_) {
        idx = pidx;
        d = 0.0;
        n_em_iter = 0;
        i1 = i1_;
        i2 = i2_;
        // pars = pars_;
        sharedSites = (int *)malloc(_sharedSites * sizeof(int));
        for (size_t i = 0; i < _sharedSites; i++) {
            sharedSites[i] = -1;
        }

        optim_jointGenoCountDist = (double *)malloc(9 * sizeof(double));
        optim_jointGenoProbDist = (double *)malloc(9 * sizeof(double));
    }
    ~pairStruct() {
        FREE(sharedSites);
        FREE(optim_jointGenoCountDist);
        FREE(optim_jointGenoProbDist);
    }

    // void print(FILE *fp, bcf_hdr_t *hdr)
    // {
    // 	fprintf(fp, "%d,%d,%d,%s,%s", pars->lut_idxToInds[idx][0], pars->lut_idxToInds[idx][1], idx, hdr->samples[pars->lut_idxToInds[idx][0]], hdr->samples[pars->lut_idxToInds[idx][1]]);
    // }

    // void print(FILE *fp)
    // {
    // 	fprintf(fp, "%d,%d,%d", pars->lut_idxToInds[idx][0], pars->lut_idxToInds[idx][1], idx);
    // }

    void sharedSites_add(size_t site_i) {
        if (snSites > _sharedSites) {
            sharedSites_expand();
        }
        sharedSites[snSites] = site_i;
        snSites++;
    }

    void sharedSites_expand() {
        size_t oldSize = _sharedSites;

        _sharedSites *= 2;
        sharedSites = (int *)realloc(sharedSites, _sharedSites * sizeof(int));
        ASSERT(NULL != sharedSites);
        for (size_t i = oldSize; i < _sharedSites; i++) {
            sharedSites[i] = -1;
        }
    }
};

/**
 * @brief distanceMatrixStruct stores the distance matrix
 *
 * @param M 			distance matrix
 * @param itemLabels 	labels of the items in the distance matrix
 * @param nInd 			number of individuals
 * @param nIndCmb		number of individual combinations
 * @param isSquared 	1 if the distance matrix is squared, 0 otherwise
 *
 */
struct distanceMatrixStruct {
    double *M = NULL;

    // idx2inds[pair_index][0] = index of the first individual in the pair
    // idx2inds[pair_index][1] = index of the second individual in the pair
    int **idx2inds = NULL;

    // inds2idx[i1][i2] = index of the pair (i1,i2) in the distance matrix
    int **inds2idx = NULL;

    char **itemLabels = NULL;

    int nInd = 0;
    int nIndCmb = 0;
    int isSquared = -1;

    distanceMatrixStruct(int nInd_, int nIndCmb_, int isSquared_, char **itemLabels_);
    ~distanceMatrixStruct();

    void print(int printMatrix, IO::outputStruct *out_dm_fs);

    void set_item_labels(char **indNames);
    // TODO print in csv format
    //  for (int i1=0; i1<metadataSt->nInd-1; i1++)
    //  {
    //  	for (int i2=i1+1; i2<metadataSt->nInd; i2++)
    //  	{
    //  		int pidx = nCk_idx(metadataSt->nInd, i1, i2);
    //  		fprintf(stdout,"%d,%d,%d,", i1,i2,pidx);
    //  		fprintf(stdout,"%.*f\n",DBL_MAXDIG10, dMS[0]->M[pidx]);
    //  	}
    //  }
};

/// @brief read distance matrix from distance matrix csv file
/// @param in_dm_fp input distance matrix file
/// @param pars paramStruct parameters
/// @return distance matrix double*
distanceMatrixStruct *distanceMatrixStruct_read(paramStruct *pars, argStruct *args);

// void print_SFS_GT(const char *TYPE, IO::outputStruct *out_sfs_fs, paramStruct *pars, int *SFS_GT3, int snSites, const char *sample1, const char *sample2);

struct metadataStruct {
    // total number of individuals in the entire dataset
    int nInd = 0;

    // number of hierarchical levels excluding the lowest level (i.e. individual)
    int nLevels = 0;

    int **nIndPerStrata = NULL;

    // Individual to Group Bitset Association
    // --------------------------------------
    //
    // ind, region, population
    // ind1, reg1, pop1
    // ind2, reg1, pop1
    // ind3, reg1, pop2
    // ind4, reg2, pop3
    //
    // association:
    // 		reg1, reg2, pop1, pop2, pop3
    // ind1,1,    0,    1,    0,    0
    // ind2,1,    0,    1,    0,    0
    // ind3,1,    0,    0,    1,    0
    // ind4,1,    1,    0,    0,    1

    // groupKeys[nBits]
    // access: groupKeys[bit] = group key for the group represented by 'bit'th bit
    // reg1 = 10000
    // reg2 = 01000
    // pop1 = 10100 // member of reg1
    // pop2 = 10010 // member of reg1
    // pop3 = 01001 // member of reg2
    uint64_t *groupKeys = NULL;

    // indKeys[nInd]
    // for ind3, which is from reg1 and pop2
    // ind3's key = pop2's key stored at ind3's index
    uint64_t *indKeys = NULL;

    // nGroups[nLevels]
    // access: nGroups[h_i] = number of groups at level h_i
    int *nGroups = NULL;

    // (lvl, g) -> lvlg_idx lut
    // lvl = hierarchical level
    // g = group index at level lvl
    // e.g. {reg1,reg2,pop1,pop2,pop3}
    // 		(0,1) -> 1 (lvl 0, group at index 1 in lvl 0)==reg2
    // 		(1,2) -> 4 (lvl 1, group at index 2 in lvl 1)==pop3
    // lvlgToIdx[lvl][g] = idx
    // idxToLvlg[idx][0] = lvl
    // idxToLvlg[idx][1] = g
    int **lvlgToIdx = NULL;
    int **idxToLvlg = NULL;

    // indNames[i] = name of individual i (char *)
    char **indNames = NULL;

    // groupNames[max_n_levels][max_n_groups_per_level][max_group_name_length]
    // access: groupNames[0][0] = "group1"
    char ***groupNames = NULL;

    // levelNames[nLevels+1] = {individual, level1, level2, ..., level_n}
    // names of levels in the hierarchy (e.g. region, population, subpopulation, individual)
    // NOTE: levelNames[0] = "individual", therefore the indexing is shifted by +1
    char **levelNames = NULL;

    // lvlStartPos[lvl] = index of the first bit in the group key corresponding to the group at level lvl
    int *lvlStartPos = NULL;

    // [nBits] total number of bits used in the construction of individual keys
    // e.g. metadata table
    // ind1, reg1, pop1, subpop1
    // ind2, reg1, pop1, subpop1
    // ind3, reg1, pop2, subpop2
    // ind4, reg2, pop2, subpop3
    // nBits = size of {reg1, reg2} + size of {pop1, pop2} + size of {subpop1, subpop2, subpop3} = 2 + 2 + 3 = 7
    int nBits = 0;

    metadataStruct(int nInd);
    ~metadataStruct();

    void print_indKeys();
    void print_groupKeys();
    void printAll();

    void resize();

    uint64_t get_indKey(int ind_i) {
        ASSERT(ind_i < nInd);
        return groupKeys[lvlgToIdx[nLevels - 1][indKeys[ind_i]]];
    }

    int get_lvlgidx(int lvl, int g) {
        return (lvl * nGroups[lvl - 1] + g);
    }

    /// @param lvl:		hierarchical level of the group to add
    /// @param g:	group index at level lvl of the group to add
    /// @param name:	name of the group to add
    void addGroup(int lvl, int g, char *name) {
        IO::vprint(2, "Found new group: %s at level %ld with index %ld", name, lvl, g);

        // add group name
        groupNames[lvl][g] = (char *)malloc(sizeof(char) * (strlen(name) + 1));
        ASSERT(strncpy(groupNames[lvl][g], name, strlen(name) + 1) != NULL);
        nGroups[lvl]++;
    }

    void setGroupKey(int bit_i, int lvl, int g, int prev_bit) {
        // fprintf(stderr, "setting the key for group %s at level %d (bit %d) with parent %d (bit %d) with key %ld (bit %d)\n", groupNames[lvl][g], lvl, bit_i, lvl-1, prev_bit, groupKeys[prev_bit], prev_bit);
        // DEVPRINT("groupName[%d][%d] = %s", lvl, g, groupNames[lvl][g]);

        // if not the first group at this level (prev_bit != -1)
        // prev_bit is the bit_i of the parent group
        // carry the parent's key at groupKeys[prev_bit] to the child group at bit_i

        if (prev_bit != -1) {
            ASSERT(groupKeys[prev_bit] > 0);
            groupKeys[bit_i] = groupKeys[prev_bit];
            // DEVPRINT("group %s at level %d group_idx_in_level %d has a parent %s at level %d group_idx_in_level %d. the parent's key is set to %ld. the group has a key of %ld before setting it", groupNames[lvl][g], lvl, g, groupNames[parent_lvl][parent_g], parent_lvl, parent_g , groupKeys[prev_bit], groupKeys[bit_i]);
        } else {
            groupKeys[bit_i] = 0;
            // DEVPRINT("group %s at level %d group_idx_in_level %d without parent. the group has a key of %ld before setting it", groupNames[lvl][g], lvl, g, groupKeys[bit_i]);
        }

        // set the bit assigned for representing this group
        // if no parent, it starts from 0 as all keys are initialized to 0 during construction
        // if parent, it starts from the parent's key as the parent's key is carried to the child
        BITSET(groupKeys[bit_i], bit_i);
        // DEVPRINT("group %s after setting the key %ld", groupNames[lvl][g], groupKeys[bit_i]);
    }

    /// @brief addLevelName - add a level name to the metadata structure
    /// @param levelName  - name of the level
    /// @param level_idx  - index of the level
    void addLevelName(const char *levelName, const int level_idx);

    // discard the bits lower than the level at interest
    // e.g.  (lvl3 \isSubsetOf lvl2 \isSubsetOf lvl1)
    //
    // number of groups at level 1: 4 == |{reg1,reg2,reg3,reg4}|
    // number of groups at level 2: 3 == |{pop1,pop2,pop3}|
    // number of groups at level 3: 5 == |{subpop1,subpop2,subpop3,subpop4,subpop5}|
    //
    // assume pop3 is a subset of reg1
    // lvlStartPos={0,4,7}
    //
    // if we are checking if ind4 belongs to pop3, which is a group at level 2(1based)
    // we discard the bits representing the groups at level 3(1based)
    //
    // pop3 key:
    // 1000 0010 0000 0000 0000 0000 0000 0000
    // ^	^  ^
    // |	|  |_ 7 == lvlStartPos for level 3(1based)
    // |    |_ 4 == lvlStartPos for level 2(1based)
    // |_ 0 == lvlStartPos for level 1(1based)
    //
    //
    // assume ind4 key: (ind4 is from subpop5, subpop5 is from pop3, pop3 is from reg1)
    // 1000 0010 0001 0000 0000 0000 0000 0000
    // ^	^  ^==============================
    //			we can discard these bits starting from the 7th bit
    //

    /// @brief indsFromGroup - check if both of the given individuals belong to a given group
    /// @param ind1
    /// @param ind2
    /// @param globGrpIdx
    /// @return
    int indsFromGroup(int ind1, int ind2, int globGrpIdx);

    // TODO
    //  int indPairFromGroup(int pair_idx, int globGrpIdx);

    /// @brief countIndsInGroup - count the number of individuals in a given group
    /// @param lvl  - hierarchical level of the group
    /// @param localGrpIdx - group index at level lvl (local to the level)
    /// @return  - number of individuals in the group
    int countIndsInGroup(int lvl, int localGrpIdx);

    /// @brief countIndsInGroup - count the number of individuals in a given group
    /// @param globIdx - group index (global, == its bit)
    /// @return  - number of individuals in the group
    int countIndsInGroup(int globIdx);

    /// @brief indFromGroup - check if a given individual belongs to a given group
    /// @param ind_i
    /// @param lvl_i
    /// @param localGrpIdx
    /// @return 1 if the individual belongs to the group, 0 otherwise
    int indFromGroup(int ind_i, int lvl_i, int localGrpIdx);

    /// @brief getNIndPerStrata - count the number of individuals in each strata
    /// and save it in the nIndPerStrata array
    void getNIndPerStrata();

    /// @brief groupFromParentGroup - check if a given group is a children of a given parent group at a given level
    /// @param plvl		- parent level
    /// @param pg		- parent group index (local)
    /// @param lvl 		- level of the children group to check
    /// @param g		- children group index (local)
    /// @return			1 if the group is a children of the parent group, 0 otherwise
    int groupFromParentGroup(int plvl, int pg, int lvl, int g);

    /// @brief countNSubgroupAtLevel - count the number of subgroups a given group has at a given level
    /// @param plvl
    /// @param pg
    /// @param lvl
    /// @return
    int countNSubgroupAtLevel(int plvl, int pg, int lvl);

    /// @brief whichLevel - get the 1-based level index of a given level name
    /// @param levelName  - name of the level
    /// @return index of the level, throw an error if the level name is not found
    int whichLevel1(const char *levelName);
};

metadataStruct *metadataStruct_get(argStruct *args, paramStruct *pars);
void metadataStruct_destroy(metadataStruct *mtd);

struct indPairThreads {
    pairStruct *pair;
    double **lngls;

    double tole = 0.0;
    int maxEmIter = 0;

    size_t nSites;

    indPairThreads(pairStruct *tPair, double **lngl, argStruct *args, paramStruct *pars) {
        pair = tPair;
        lngls = lngl;
        // args = args_;
        // pars = pars_;
        tole = args->tole;
        maxEmIter = args->maxEmIter;
        nSites = pars->nSites;
        // ASSERT(pars->nSites > 0);
        // tes = pars->nSites;
    }
};

// struct bootstrapThreads {
//     amovaStruct **amova;
//     double **lngls;
//     size_t nSites;

//     bootstrapThreads();
//     ~bootstrapThreads();
// };

#endif  // __DATA_STRUCTS__