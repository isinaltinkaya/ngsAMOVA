#ifndef __DATA_STRUCTS__
#define __DATA_STRUCTS__

#include "argStruct.h"
#include "mathUtils.h"
#include "paramStruct.h"
#include "vcfReader.h"





/* FORWARD DECLARATIONS ----------------------------------------------------- */

namespace IO {
    void validateString(const char* str);
    int verbose(const int verbose_threshold);
    void vprint(const char* format, ...);
    void vprint(const int verbose_threshold, const char* format, ...);
    void vprint(FILE* fp, const int verbose_threshold, const char* format, ...);
}

struct blobStruct;

typedef struct vcfData vcfData;
typedef struct indPairThreads indPairThreads;

typedef struct intArray intArray;

/// ----------------------------------------------------------------------- ///

struct intArray {

    /// @var d - array of int data
    int* d;

    /// @var size - size of d
    size_t size;

    /// @var len - length of the currently used part of d (
    /// @note 
    /// len   = position right after the last element in use in d
    ///         the next element will be added at position len
    /// len-1 = gives the index of the last element in d
    /// if intArray size is preallocated, len may be < size-1
    /// if len == 0 the data is empty
    /// if len == size the data is full
    /// len can never be > size
    size_t len;

    /// @brief add - add an int
    /// @param i int to add
    /// @return size_t index of the newly added int
    size_t add(const int i) {

        // note: can never happen since during initialization d is mallocated with size 1 and size is set to 1
        // if (0 == this->size) {
        //     this->size = 1;
        //     this->d = (int*)malloc(sizeof(int));
        //     ASSERT(this->d != NULL);
        // }

        if (this->len == this->size) {
            this->size++;
            int* tmp = NULL;
            tmp = (int*)realloc(this->d, this->size * sizeof(int));
            ASSERT(tmp != NULL);
            this->d = tmp;
        } else if (this->len > this->size) {
            NEVER;
        }


        this->d[this->len] = i;
        this->len++;
        return (this->len - 1);
    }

    bool is_full(void) {
        return(this->len == this->size);
    }

    /// @brief find - get index of the given value in the intArray
    /// @param val - int value to search for
    /// @param val_idx - pointer to the variable to store the index of the given int (if found)
    /// @return bool true if the given int is found, false otherwise
    bool find(const int val, size_t* val_idx) {
        for (size_t i = 0;i < this->len;++i) {
            if (val == this->d[i]) {
                *val_idx = i;
                return(true);
            }
        }
        return(false);
    }

    /// @brief contains - check if the given value is in the intArray
    /// @param val - int value to search for
    /// @return bool true if the given int is found, false otherwise
    bool contains(const int val) {
        for (size_t i = 0;i < this->len;++i) {
            if (val == this->d[i]) {
                return(true);
            }
        }
        return(false);
    }


    /// @brief get - get the value at the given index
    /// @param idx index of the value to get
    /// @return int value at the given index
    int get(const size_t idx) {
        ASSERT(idx < this->size);
        ASSERT(idx <= this->len);
        return (this->d[idx]);
    }

    /// @brief clear - clear the intArray
    /// @note does not free memory
    void clear(void) {
        this->len = 0;
    }

};

inline intArray* intArray_alloc(const size_t size) {
    ASSERT(size > 0);
    intArray* ret = (intArray*)malloc(sizeof(intArray));
    ASSERT(ret != NULL);
    ret->d = NULL;
    ret->d = (int*)malloc(size * sizeof(int));
    ASSERT(ret->d != NULL);
    ret->size = size;
    ret->len = 0;
    return(ret);
}

inline intArray* intArray_init(void) {
    intArray* ret = (intArray*)malloc(sizeof(intArray));
    ASSERT(ret != NULL);
    ret->d = NULL;
    ret->d = (int*)malloc(sizeof(int));
    ASSERT(ret->d != NULL);
    ret->size = 1;
    ret->len = 0;
    return(ret);
}

inline intArray* intArray_init(const int founder_val) {
    intArray* ret = (intArray*)malloc(sizeof(intArray));
    ASSERT(ret != NULL);
    ret->d = NULL;
    ret->d = (int*)malloc(sizeof(int));
    ASSERT(ret->d != NULL);
    ret->d[0] = founder_val;
    ret->size = 1;
    ret->len = 1;
    return(ret);
}

inline void intArray_destroy(intArray* ia) {
    FREE(ia->d);
    FREE(ia);
}

typedef struct size_tArray size_tArray;
struct size_tArray {

    /// @var d - array of 
    size_t* d;

    /// @var size - size of d
    size_t size;

    /// @var len - position right after the last element in use in d
    /// @note 
    /// len-1 gives the index of the last element in d
    /// if size_tArray size is preallocated, len may be < size-1
    size_t len;

    /// @brief add - add a size_t
    /// @param i size_t to add
    /// @return size_t index of the newly added size_t
    size_t add(const size_t i) {

        // if (0 == this->size) {
        //     this->size = 1;
        //     this->d = (size_t*)malloc(sizeof(size_t));
        //     ASSERT(this->d != NULL);
        // }

        if (this->len == this->size) {
            this->size++;
            size_t* tmp = NULL;
            tmp = (size_t*)realloc(this->d, this->size * sizeof(size_t));
            ASSERT(tmp != NULL);
            this->d = tmp;
            tmp = NULL;

        } else if (this->len > this->size) {
            NEVER;
        }


        this->d[this->len] = i;
        this->len++;
        return (this->len - 1);
    }


    /// @brief set_val - set the value at the given index
    /// @param val 
    /// @param idx 
    /// @note use with caution, does not check if the index is already in use, and does not track if any of the previous vals are set/unset
    void set_val(const size_t val, const size_t idx) {
        ASSERT(idx < this->size);
        ASSERT(idx >= this->len);
        this->d[idx] = val;
        if (len == idx) {
            this->len++;
        }
        return;
    }




    /// @brief find - get index of the given value in the size_tArray
    /// @param val - size_t value to search for
    /// @param val_idx - pointer to the variable to store the index of the given size_t (if found)
    /// @return bool true if the given size_t is found, false otherwise
    bool find(const size_t val, size_t* val_idx) {
        for (size_t i = 0;i < this->len;++i) {
            if (val == this->d[i]) {
                *val_idx = i;
                return(true);
            }
        }
        return(false);
    }

    /// @brief contains - check if the given value is in the size_tArray
    /// @param val - size_t value to search for
    /// @return bool true if the given size_t is found, false otherwise
    bool contains(const size_t val) {
        for (size_t i = 0;i < this->len;++i) {
            if (val == this->d[i]) {
                return(true);
            }
        }
        return(false);
    }

    /// @brief get - get the value at the given index
    /// @param idx index of the value to get
    /// @return size_t value at the given index
    size_t get(const size_t idx) {
        ASSERT(idx < this->size);
        ASSERT(idx <= this->len);
        return (this->d[idx]);
    }

    /// @brief clear - clear the size_tArray
    /// @note does not free memory
    void clear(void) {
        this->len = 0;
    }

};

inline size_tArray* size_tArray_alloc(const size_t size) {
    ASSERT(size > 0);
    size_tArray* ret = (size_tArray*)malloc(sizeof(size_tArray));
    ASSERT(ret != NULL);
    ret->d = NULL;
    ret->d = (size_t*)malloc(size * sizeof(size_t));
    ASSERT(ret->d != NULL);
    ret->size = size;
    ret->len = 0;
    return(ret);
}

inline size_tArray* size_tArray_init(void) {
    size_tArray* ret = (size_tArray*)malloc(sizeof(size_tArray));
    ASSERT(ret != NULL);
    ret->d = NULL;
    ret->d = (size_t*)malloc(sizeof(size_t));
    ASSERT(ret->d != NULL);
    ret->size = 1;
    ret->len = 0;
    return(ret);
}

/// @brief size_tArray_init - allocate memory for a size_tArray and initialize it with a founder value
inline size_tArray* size_tArray_init(const size_t founder_val) {
    size_tArray* ret = (size_tArray*)malloc(sizeof(size_tArray));
    ret->d = (size_t*)malloc(sizeof(size_t));
    ASSERT(ret->d != NULL);
    ret->size = 1;
    ret->d[0] = founder_val;
    ret->len = 1;
    return(ret);
}

inline void size_tArray_destroy(size_tArray* sa) {
    FREE(sa->d);
    FREE(sa);
}

typedef struct strArray strArray;
struct strArray {

    /// @var d - array of char* data
    char** d;

    /// @var size - size of d
    size_t size;

    /// @var len - position right after the last element in use in d
    /// @note
    /// len-1 gives the index of the last element in d
    /// if strArray size is preallocated, len may be < size-1
    size_t len;

    /// @brief add - add a string (char*)
    /// @param str string to add
    /// @return size_t index of the newly added string (char*)
    size_t add(const char* str) {

        ASSERT(str != NULL);

        if (this->len == this->size) {
            this->size++;
            char** tmp = NULL;
            DEVASSERT(NULL != this->d);
            tmp = (char**)realloc(this->d, this->size * sizeof(char*));
            ASSERT(tmp != NULL);
            this->d = tmp;
            this->d[this->len] = NULL;
        } else if (this->len > this->size) {
            NEVER;
        }

        this->d[this->len] = strdup(str);
        ASSERT(this->d[this->len] != NULL);
        this->len++;
        return (this->len - 1);
    }

    void add_to(const char* val, const size_t idx) {
        ASSERT(idx < this->size);
        ASSERT(this->size > 0);
        ASSERT(NULL != this->d);
        ASSERT(NULL == this->d[idx]);
        this->d[idx] = strdup(val);
        ASSERT(this->d[idx] != NULL);
        if (len == idx) {
            this->len++;
        } else if (len < idx) {
            this->len = idx + 1;
        }
        return;
    }

    /// @brief find - get index of the given value in the strArray
    /// @param val - str (char*) value to search for
    /// @param val_idx - pointer to the variable to store the index of the given str (if found)
    /// @return bool true if the given str is found, false otherwise
    bool find(const char* val, size_t* val_idx) {
        ASSERT(val != NULL);
        for (size_t i = 0;i < this->len;++i) {
            // can be NULL if a specific value is set using add_to which skipped 1 or more indices
            if (NULL != this->d[i]) {
                if (0 == strcmp(val, this->d[i])) {
                    *val_idx = i;
                    return(true);
                }
            }
        }
        return(false);
    }

    /// @brief find_from - get index of the given value in the strArray starting from the given index
    /// @param val - str (char*) value to search for
    /// @param val_idx - pointer to the variable to store the index of the given str (if found)
    /// @param start_idx - index to start the search from
    /// @return bool true if the given str is found, false otherwise
    bool find_from(const char* val, size_t* val_idx, const size_t start_idx) {
        ASSERT(val != NULL);
        for (size_t i = start_idx;i < this->len;++i) {
            // can be NULL if a specific value is set using add_to which skipped 1 or more indices
            if (NULL != this->d[i]) {
                if (0 == strcmp(val, this->d[i])) {
                    *val_idx = i;
                    return(true);
                }
            }
        }
        return(false);
    }

    /// @brief contains - check if the given value is in the strArray
    /// @param val - str (char*) value to search for
    /// @return bool true if the given str is found, false otherwise
    bool contains(const char* val) {
        ASSERT(val != NULL);
        DEVASSERT(this->d != NULL);
        for (size_t i = 0;i < this->len;++i) {
            // can be NULL if a specific value is set using add_to which skipped 1 or more indices
            if (NULL != this->d[i]) {
                if (0 == strcmp(val, this->d[i])) {
                    return(true);
                }
            }
        }
        return(false);
    }

    /// @brief get - get the value at the given index
    /// @param idx index of the value to get
    /// @return char* value at the given index
    char* get(const size_t idx) {
        ASSERT(idx < this->size);
        if (idx < this->len) {
            return (this->d[idx]);
        }
        return(NULL);
    }

    void clear(void) {
        for (size_t i = 0;i < this->len;++i) {
            FREE(this->d[i]);
        }
        this->len = 0;
    }

};

inline strArray* strArray_init(void) {
    strArray* ret = (strArray*)malloc(sizeof(strArray));
    ASSERT(ret != NULL);
    ret->d = NULL;
    ret->size = 1;
    ret->len = 0;
    ret->d = (char**)malloc(sizeof(char*));
    ASSERT(ret->d != NULL);
    ret->d[0] = NULL;
    return(ret);
}

inline strArray* strArray_init(const char* founder_val) {
    ASSERT(founder_val != NULL);

    strArray* ret = (strArray*)malloc(sizeof(strArray));
    ASSERT(ret != NULL);
    ret->d = NULL;
    ret->size = 1;
    ret->len = 1;

    ret->d = (char**)malloc(sizeof(char*));
    ASSERT(ret->d != NULL);
    ret->d[0] = strdup(founder_val);
    ASSERT(ret->d[0] != NULL);
    return(ret);
}

/// @brief strArray_alloc - allocate memory for a strArray
/// @param size size of the strArray to allocate
/// @return strArray* pointer to the allocated strArray
inline strArray* strArray_alloc(const size_t size) {
    ASSERT(size > 0);

    strArray* ret = (strArray*)malloc(sizeof(strArray));
    ASSERT(ret != NULL);

    ret->d = NULL;
    ret->d = (char**)malloc(size * sizeof(char*));
    ASSERT(ret->d != NULL);
    for (size_t i = 0;i < size;++i) {
        ret->d[i] = NULL;
    }
    ret->len = 0;
    ret->size = size;
    return(ret);
}

inline void strArray_destroy(strArray* sa) {
    for (size_t i = 0;i < sa->size;++i) {
        if (sa->d[i] != NULL) {
            FREE(sa->d[i]);
        }
    }
    FREE(sa->d);
    FREE(sa);
}





/* -------------------------------------------------------------------------- */


/// @brief read_formula_str - read the formula string into a strArray of formula tokens
/// @param formula formula string
/// @return pointer to newly created strArray
/// @example formula = 'Individual ~ Region/Population/Subpopulation'
/// formulaTokens->len = 4
/// formulaTokens->d[0] = "Region" (level 1, lvlidx 0)
/// formulaTokens->d[1] = "Population" (level 2, lvlidx 1)
/// formulaTokens->d[2] = "Subpopulation" (level 3, lvlidx 2)
/// formulaTokens->d[3] = "Individual" (level 4, lvlidx 3)
strArray* read_formula_str(const char* formula);

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
struct metadataStruct {

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
    /// @note free iff pars->metadata==NULL && PROGRAM_HAS_INPUT_VCF ; else it is ptr to pars->metadata->indNames
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

/// @param n_levels number of hierarchical levels in the formula (i.e. formulaTokens->len)
inline metadataStruct* metadataStruct_init(const int in_nLevels) {

    metadataStruct* ret = (metadataStruct*)malloc(sizeof(metadataStruct));
    ASSERT(ret != NULL);

    ret->nLevels = in_nLevels;
    ret->nGroups = 0;
    ret->nInd = 0;

    ret->levelNames = NULL;
    ret->levelNames = strArray_alloc(in_nLevels);

    ret->indNames = NULL;
    ret->indNames = strArray_init();

    ret->groupNames = NULL;
    ret->groupNames = strArray_init();

    size_t nLevelsExclInd = in_nLevels - 1;

    ret->level2groupIndices = NULL;
    ret->level2groupIndices = (size_tArray**)malloc(nLevelsExclInd * sizeof(size_tArray*));
    ASSERT(ret->level2groupIndices != NULL);
    for (size_t i = 0;i < (size_t)nLevelsExclInd;++i) {
        ret->level2groupIndices[i] = NULL;
        ret->level2groupIndices[i] = size_tArray_alloc(1);
    }

    ret->group2levelIndices = NULL;
    ret->group2levelIndices = size_tArray_alloc(1);



    ret->group2indIndices = NULL;
    ret->group2indIndices = (size_tArray**)malloc(sizeof(size_tArray*));
    ASSERT(ret->group2indIndices != NULL);
    ret->group2indIndices[0] = NULL;
    ret->group2indIndices[0] = size_tArray_alloc(1);

    ret->group2pairIndices = NULL;
    ret->group2pairIndices = (size_tArray**)malloc(sizeof(size_tArray*));
    ASSERT(ret->group2pairIndices != NULL);
    ret->group2pairIndices[0] = NULL;
    ret->group2pairIndices[0] = size_tArray_alloc(1);


    ret->group2subgroupIndices = NULL;
    ret->group2subgroupIndices = (size_tArray**)malloc(sizeof(size_tArray*));
    ASSERT(ret->group2subgroupIndices != NULL);
    ret->group2subgroupIndices[0] = NULL;
    ret->group2subgroupIndices[0] = size_tArray_alloc(1);

    return(ret);

}


/// @brief read metadata from metadata file into metadataStruct
metadataStruct* metadataStruct_read(paramStruct* pars);



inline void metadataStruct_destroy(metadataStruct* mtd) {

    strArray_destroy(mtd->levelNames);
    strArray_destroy(mtd->indNames);
    strArray_destroy(mtd->groupNames);


    for (size_t i = 0;i < (size_t)mtd->nLevels - 1;++i) {
        size_tArray_destroy(mtd->level2groupIndices[i]);
    }
    FREE(mtd->level2groupIndices);

    for (size_t i = 0;i < (size_t)mtd->nGroups;++i) {

        size_tArray_destroy(mtd->group2indIndices[i]);

        size_tArray_destroy(mtd->group2pairIndices[i]);
        if (NULL != mtd->group2subgroupIndices[i]) {
            // is null if the group is in the last level
            // thus effective size in use for this array is actually nGroups-nGroupsAtLastLevel
            size_tArray_destroy(mtd->group2subgroupIndices[i]);
        }
    }

    size_tArray_destroy(mtd->group2levelIndices);

    FREE(mtd->group2indIndices);
    FREE(mtd->group2pairIndices);
    FREE(mtd->group2subgroupIndices);

    FREE(mtd);
}


struct indPairThreads {

    paramStruct* pars = NULL;
    vcfData* vcfd = NULL;

    int pidx = -1;

    indPairThreads(vcfData* vcfd_, paramStruct* pars_, const int pairIndex) {
        this->pars = pars_;
        this->vcfd = vcfd_;
        pidx = pairIndex;
    }
};

inline void trimSpaces(char* str) {
    char* ptr = str;
    char* end = NULL;

    // skip leading spaces and find the start of the non-space chars
    while (*ptr && isspace((unsigned char)*ptr)) {
        ++ptr;
    }

    if (*ptr == '\0') {
        // it is all spaces!?
        ERROR("Found empty string.");
    }

    // find the end of the str and set end to the position of last non-space char
    for (char* tmp = ptr; *tmp; ++tmp) {
        if (!isspace((unsigned char)*tmp)) {
            end = tmp;
        }
    }

    if (ptr != str) {
        // we have leading spaces
        // shift the non-space chars to the beginning
        char* dest = str;
        while (ptr <= end) {
            *dest++ = *ptr++;
        }
        *dest = '\0'; // terminate the shifted str
    } else {
        // no leading spaces
        // terminate the str at the last non-space char
        *(end + 1) = '\0';
    }
}

typedef struct jgtmat_t jgtmat_t;
struct jgtmat_t {

    // typically: nIndCmb (one matrix per ind pair)
    // or nIndCmb*(args->nBootstraps+1) (for each bootstrap replicate+the orig)
    size_t n;      // number of matrices

    // - 9 (3x3) for MM Mm mm
    // - 100 (10x10) for all genotype combinations
    size_t size;   // size of each matrix

    // probabilities
    // pm[n][size]
    double** pm;

    // counts
    // m[n][size]
    uint64_t** m;

};

inline jgtmat_t* jgtmat_init(const size_t nmat) {
    jgtmat_t* ret = NULL;
    ret = (jgtmat_t*)malloc(sizeof(jgtmat_t));
    ASSERT(ret != NULL);

    ret->n = nmat;
    ret->pm = NULL;
    ret->m = NULL;

    if (ARG_DOJGTM_3GT == args->doJGTM) {
        ret->size = 9;
    } else if (ARG_DOJGTM_10GT == args->doJGTM) {
        ret->size = 100;
    } else {
        NEVER;
    }

    if (args->doEM) {
        ret->pm = (double**)malloc(ret->n * sizeof(double*));
        ASSERT(ret->pm != NULL);
        for (size_t i = 0;i < ret->n;++i) {
            ret->pm[i] = NULL;
            ret->pm[i] = (double*)malloc(ret->size * sizeof(double));
            ASSERT(ret->pm[i] != NULL);
            if (9 == ret->size) {
                // use flat prior
                ret->pm[i][0] = FRAC_1_9;
                ret->pm[i][1] = FRAC_1_9;
                ret->pm[i][2] = FRAC_1_9;
                ret->pm[i][3] = FRAC_1_9;
                ret->pm[i][4] = FRAC_1_9;
                ret->pm[i][5] = FRAC_1_9;
                ret->pm[i][6] = FRAC_1_9;
                ret->pm[i][7] = FRAC_1_9;
                ret->pm[i][8] = FRAC_1_9;
            } else if (100 == ret->size) {
                // use flat prior
                for (size_t j = 0;j < ret->size;++j) {
                    ret->pm[i][j] = 0.01;
                }
            }
        }
    } else {
        ret->m = (uint64_t**)malloc(ret->n * sizeof(uint64_t*));
        ASSERT(ret->m != NULL);
        for (size_t i = 0;i < ret->n;++i) {
            ret->m[i] = NULL;
            ret->m[i] = (uint64_t*)malloc(ret->size * sizeof(uint64_t));
            ASSERT(ret->m[i] != NULL);
            for (size_t j = 0;j < ret->size;++j) {
                ret->m[i][j] = 0;
            }

        }
    }

    return(ret);
}

inline void jgtmat_destroy(jgtmat_t* d) {
    if (d->m != NULL) {
        for (size_t i = 0;i < d->n;++i) {
            FREE(d->m[i]);
        }
        FREE(d->m);
    }
    if (d->pm != NULL) {
        for (size_t i = 0;i < d->n;++i) {
            FREE(d->pm[i]);
        }
        FREE(d->pm);
    }
    FREE(d);
    return;
}

void jgtmat_get_srcgl(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd);
void jgtmat_get_srcgt(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd);



/// ----------------------------------------------------------------------- ///


#endif  // __DATA_STRUCTS__