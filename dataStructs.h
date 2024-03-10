#ifndef __DATA_STRUCTS__
#define __DATA_STRUCTS__

#include "argStruct.h"
#include "mathUtils.h"
#include "paramStruct.h"
#include "vcfReader.h"


/// ----------------------------------------------------------------------- ///



/// TODO !! MATCH THE INDNAMES IN THE VCF FILE WITH THE INDNAMES IN THE METADATA FILE


/* FORWARD DECLARATIONS ----------------------------------------------------- */

namespace IO {
    void validateString(const char* str);
    int verbose(const int verbose_threshold);
    void vprint(const char* format, ...);
    void vprint(const int verbose_threshold, const char* format, ...);
    void vprint(FILE* fp, const int verbose_threshold, const char* format, ...);
    void vvprint(FILE* fp, const int verbose_threshold, const char* format, ...);
}

struct blobStruct;

typedef struct lnglStruct lnglStruct;
typedef struct vcfData vcfData;
typedef struct indPairThreads indPairThreads;

typedef struct intArray intArray;

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

struct lnglStruct {

    /*
     * lngl[nSites][nInd*nGT*double]
     * genotype likelihoods in natural log
     * for each individual at each site
     *
     *
     */
     // d[size1][size2]
    double** d = NULL;

    // size1 follows pars->nSites_arrays_size
    size_t size1 = -1;

    size_t size2 = -1;

    lnglStruct(const int nGt, const int nInd);
    ~lnglStruct();

    void expand();

};

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

    // \def names->d[nInds]
    // names->d[i] = name of individual i (char *)
    strArray* names;

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

    ret->names = NULL;
    ret->names = strArray_init();

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
    strArray_destroy(mtd->names);
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

/// ----------------------------------------------------------------------- ///
/// DISTANCE MATRIX

/// @brief DMAT_TYPE_* - type of the distance matrix
///
///  ___  ___  ___  (exclude=0,include=1)
///   .    .  [0|1] -> diagonal
///   .  [0|1]  .   -> lower triangular
/// [0|1]  .    .   -> upper triangular
///
#define DMAT_TYPE_INCLUDE_DIAG (1<<0)
#define DMAT_TYPE_INCLUDE_LOWER_TRI (1<<1)
#define DMAT_TYPE_INCLUDE_UPPER_TRI (1<<2)
/// LTED: lower triangular, excluding the diagonal
/// (0<<2)|(1<<1)|(0<<0)
#define DMAT_TYPE_LTED 2
/// LTID: lower triangular, including the diagonal
/// (0<<2)|(1<<1)|(1<<0)
#define DMAT_TYPE_LTID 3
/// UTED: upper triangular, excluding the diagonal
/// (1<<2)|(0<<1)|(0<<0)
#define DMAT_TYPE_UTED 4
/// UTID: upper triangular, including the diagonal
/// (1<<2)|(0<<1)|(1<<0)
#define DMAT_TYPE_UTID 5
/// FULL: full matrix
/// (1<<2)|(1<<1)|(1<<0)
#define DMAT_TYPE_FULL 7


/// @brief DMAT_TRANSFORM_* - transformation applied to the distances in the matrix
///
/// NONE: not transformed
#define DMAT_TRANSFORM_NONE 0
/// SQUARE: val^2
#define DMAT_TRANSFORM_SQUARE 1
/// SQRT: sqrt(val) //TODO
// #define DMAT_TRANSFORM_SQRT 2 
/// ABS: absolute value (|val|) //TODO
// #define DMAT_TRANSFORM_ABS 3
/// PSEUDO_EUCLIDEAN: // TODO
// #define DMAT_TRANSFORM_PSEUDO_EUCLIDEAN 4


/// @brief DMAT_METHOD_* - method used in distance calculation (i.e. distance metric)
#define DMAT_METHOD_DIJ 0
#define DMAT_METHOD_SIJ 1
#define DMAT_METHOD_FIJ 2
#define DMAT_METHOD_IBS0 3
#define DMAT_METHOD_IBS1 4
#define DMAT_METHOD_IBS2 5
#define DMAT_METHOD_KIN 6
#define DMAT_METHOD_R0 7
#define DMAT_METHOD_R1 8

/// @brief DMAT_NAMES_SRC_* - source of the names in the distance matrix
/// NONE: no names (names=NULL)
#define DMAT_NAMES_SRC_NONE 0
/// IN_DM_FILE: names is allocated and read from the distance matrix input file
/// <requires cleaning>
#define DMAT_NAMES_SRC_IN_DM_FILE 1
/// IN_VCF_FILE: names is pointer to the strArray* names in pars which was read from the VCF file
/// <no cleaning>
#define DMAT_NAMES_SRC_IN_VCF_PARS_PTR 2
/// IN_METADATA_FILE: names is pointer to the strArray* names in metadataStruct
/// <no cleaning>
/// NOTE: currently not used
#define DMAT_NAMES_SRC_METADATA_NAMES_PTR 3

typedef struct dmat_t dmat_t;

/// @struct dmat_t - distance matrix struct 
/// @brief  struct for n distance matrices 
struct dmat_t {

    /// @var type - type of the distance matrix 
    /// @details bitset
    uint8_t type;

    /// @var transform - transformation applied to the distances in the matrix
    uint32_t transform;

    /// @var method - method used to calculate the distances in the matrix
    /// @details
    /// 0: Dij
    /// 1: Sij
    /// 2: Fij
    /// 3: IBS0
    /// 4: IBS1
    /// 5: IBS2
    /// 6: Kin
    /// 7: R0
    /// 8: R1
    uint32_t method;

    /// @var size - size of a each matrix matrix[n]
    /// @details typically nIndCmb
    size_t size;

    /// @var names - array of names of the individuals in the distance matrix
    /// @details
    /// - if names_src == DMAT_NAMES_SRC_IN_DM_FILE, names is allocated and read from the distance matrix input file
    /// - if names_src == DMAT_NAMES_SRC_IN_VCF_PARS_PTR, names is pointer to the strArray* names in pars which was read from the VCF file
    /// - if names_src == DMAT_NAMES_SRC_METADATA_NAMES_PTR, names is pointer to the strArray* names in metadataStruct
    /// - if names_src == DMAT_NAMES_SRC_NONE, names is NULL
    strArray* names;
    uint8_t names_src;

    /// @var n - number of distance matrices
    size_t n;

    /// @var matrix - array of n 1D distance matrices of size 'size'
    // matrix[n][size]
    double** matrix;

};


inline dmat_t* dmat_init(const size_t nInd, const uint8_t type, const uint32_t method, const uint32_t transform, strArray* names, const uint8_t names_src) {

    dmat_t* ret = NULL;
    ret = (dmat_t*)malloc(sizeof(dmat_t));
    ASSERT(ret != NULL);

    switch (type) {
    case DMAT_TYPE_LTED:
    case DMAT_TYPE_UTED:
        ret->size = (nInd * (nInd - 1)) / 2;
        break;
    case DMAT_TYPE_LTID:
    case DMAT_TYPE_UTID:
        ret->size = (nInd * (nInd + 1)) / 2;
        break;
    case DMAT_TYPE_FULL:
        ret->size = nInd * nInd;
        break;
    default:
        NEVER;
    }

    ret->type = type;
    ret->transform = transform;
    ret->method = method;

    ret->n = (args->nBootstraps > 0) ? (1 + args->nBootstraps) : 1;

    ret->names_src = names_src;
    if (ret->names_src == DMAT_NAMES_SRC_IN_VCF_PARS_PTR) {
        ret->names = names;
    } else if (ret->names_src == DMAT_NAMES_SRC_METADATA_NAMES_PTR) {
        ret->names = names;
    } else {
        NEVER;
    }

    ret->matrix = NULL;
    ret->matrix = (double**)malloc(ret->n * sizeof(dmat_t*));
    ASSERT(ret->matrix != NULL);

    for (size_t i = 0; i < ret->n;++i) {
        ret->matrix[i] = NULL;
        ret->matrix[i] = (double*)malloc(ret->size * sizeof(double));
        ASSERT(ret->matrix[i] != NULL);
        for (size_t j = 0;j < ret->size;++j) {
            ret->matrix[i][j] = 0.0;
        }
    }
    return(ret);
}



inline void dmat_destroy(dmat_t* d) {
    for (size_t i = 0; i < d->n;++i) {
        FREE(d->matrix[i]);
    }
    FREE(d->matrix);

    if (d->names_src == DMAT_NAMES_SRC_IN_DM_FILE) {
        strArray_destroy(d->names);
    } else if (d->names_src == DMAT_NAMES_SRC_IN_VCF_PARS_PTR) {
        d->names = NULL;
    } else if (d->names_src == DMAT_NAMES_SRC_METADATA_NAMES_PTR) {
        d->names = NULL;
    } else if (d->names_src == DMAT_NAMES_SRC_NONE) {
        NEVER;
    } else {
        NEVER;
    }

    FREE(d);
    return;
}

dmat_t* dmat_read(const char* in_dm_fn, const uint32_t required_transform, metadataStruct* metadata);

void dmat_print(dmat_t* dmat);




typedef struct jgtmat_t jgtmat_t;
struct jgtmat_t {

    // typically: nIndCmb (one matrix per ind pair)
    size_t n;      // number of matrices

    // - 9 (3x3) for MM Mm mm
    // - 100 (10x10) for all genotype combinations
    size_t size;   // size of each matrix

    // probabilities
    // pm[n][size]
    double** pm;

    // expected counts (from doEM with gl source)
    // em[n][size]
    double** em;

    // counts
    // m[n][size]
    int** m;

};

inline jgtmat_t* jgtmat_init(const size_t nIndCmb) {
    jgtmat_t* ret = NULL;
    ret = (jgtmat_t*)malloc(sizeof(jgtmat_t));
    ASSERT(ret != NULL);

    ret->n = nIndCmb;
    ret->pm = NULL;
    ret->em = NULL;
    ret->m = NULL;


    if (args->doJGTM == 1) {
        ret->size = 9;
    } else if (args->doJGTM == 2) {
        ret->size = 100;
    }

    ret->pm = (double**)malloc(ret->n * sizeof(double*));
    ASSERT(ret->pm != NULL);
    for (size_t i = 0;i < ret->n;++i) {
        ret->pm[i] = NULL;
        ret->pm[i] = (double*)malloc(ret->size * sizeof(double));
        ASSERT(ret->pm[i] != NULL);
        for (size_t ii = 0;ii < ret->size;++ii) {
            ret->pm[i][ii] = 0.0;
        }
    }
    if (args->doEM) {
        // TODO are counts really necessary??
        ret->em = (double**)malloc(ret->n * sizeof(double*));
        ASSERT(ret->em != NULL);
        for (size_t i = 0;i < ret->n;++i) {
            ret->em[i] = NULL;
            ret->em[i] = (double*)malloc(ret->size * sizeof(double));
            ASSERT(ret->em[i] != NULL);
            if (9 == ret->size) {
                // use flat prior
                ret->em[i][0] = FRAC_1_9;
                ret->em[i][1] = FRAC_1_9;
                ret->em[i][2] = FRAC_1_9;
                ret->em[i][3] = FRAC_1_9;
                ret->em[i][4] = FRAC_1_9;
                ret->em[i][5] = FRAC_1_9;
                ret->em[i][6] = FRAC_1_9;
                ret->em[i][7] = FRAC_1_9;
                ret->em[i][8] = FRAC_1_9;
            } else if (100 == ret->size) {
                // use flat prior
                for (size_t j = 0;j < ret->size;++j) {
                    ret->em[i][j] = 0.01;
                }
            }
        }
    } else {
        ret->m = (int**)malloc(ret->n * sizeof(int*));
        ASSERT(ret->m != NULL);
        for (size_t i = 0;i < ret->n;++i) {
            ret->m[i] = NULL;
            ret->m[i] = (int*)malloc(ret->size * sizeof(int));
            ASSERT(ret->m[i] != NULL);
            for (size_t ii = 0;ii < ret->size;++ii) {
                ret->m[i][ii] = 0;
            }
        }
    }

    return(ret);
}

inline void jgtmat_destroy(jgtmat_t* d) {
    if (d->em != NULL) {
        for (size_t i = 0;i < d->n;++i) {
            FREE(d->em[i]);
        }
        FREE(d->em);
    }
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

void jgtmat_get_srcgl(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd, blobStruct* blob);
void jgtmat_get_srcgt(jgtmat_t* jgtm, paramStruct* pars, vcfData* vcfd, blobStruct* blob);



/// ----------------------------------------------------------------------- ///


#endif  // __DATA_STRUCTS__