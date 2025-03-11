/**
 * @file    dmat.h
 * @brief   header file for dmat.cpp
 * @details contains distance matrix data structure (dmat_t) and associated functions
 */
#ifndef __DMAT_H__
#define __DMAT_H__

 /* INCLUDES ----------------------------------------------------------------- */
#include "shared.h"
/* END-OF-INCLUDES ---------------------------------------------------------- */

/* FORWARD-DECLARATIONS ----------------------------------------------------- */
typedef struct metadata_t metadata_t;
typedef struct strArray strArray;
typedef struct jgtmat_t jgtmat_t;
typedef struct dmat_t dmat_t;
typedef struct outfile_t outfile_t;
/* END-OF-FORWARD-DECLARATIONS ---------------------------------------------- */

/* MACROS ------------------------------------------------------------------- */
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



#define DMAT_INTPLUS_TRANSFORM_NONE (0<<0)
#define DMAT_INTPLUS_TRANSFORM_SQUARE (1<<0)
#define DMAT_INTPLUS_TRANSFORM_CAILLIEZ (1<<1)
#define DMAT_INTPLUS_TRANSFORM_LINGOES (1<<2)
#define DMAT_INTPLUS_TRANSFORM_QUASIEUCLID (1<<3)
#define DMAT_INTPLUS_TRANSFORM_MAX (1<<1) // lingues and quasieuclid not implemented yet

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
#define DMAT_METHOD_DXY 9

const char* get_dmat_method_str(const int method);

/// @brief DMAT_NAMES_SRC_* - source of the names in the distance matrix
/// NONE: no names (names=NULL)
#define DMAT_NAMES_SRC_NONE 0
/// IN_DM_FILE: names is allocated and read from the distance matrix input file, implies input is distance matrix AND metadata is NOT provided
/// <requires cleaning>
#define DMAT_NAMES_SRC_IN_DM_FILE 1
/// IN_VCF_FILE: names is pointer to the strArray* names in pars which was read from the VCF file, implies no metadata is provided
/// <no cleaning>
#define DMAT_NAMES_SRC_IN_VCF_PARS_PTR 2
/// IN_METADATA_FILE: names is pointer to the strArray* names in metadata, implies input is distance matrix and metadata file is provided
/// <no cleaning>
#define DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR 3
/// PRIVATE: names is allocated and used internally in the program
#define DMAT_NAMES_SRC_PRIVATE 4
/* END-OF-MACROS ------------------------------------------------------------ */


/* TYPEDEF-STRUCTS ---------------------------------------------------------- */
/// @struct dmat_t - distance matrix struct 
/// @brief  struct for n distance matrices 
struct dmat_t {

    /// @var n - number of distance matrices
    /// @details 1 if no bootstrapping
    ///            args->nBootstraps + 1 (original) if bootstrapping enabled
    size_t n;

    /// @var size - size of a each matrix matrix[n]
    /// @details typically nIndCmb
    size_t size;

    /// @var type - type of the distance matrix 
    /// @details bitset
    uint8_t type;

    /// @var transform - transformation applied to the distances in the matrix
    uint8_t transform;

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
    /// 9: Dxy
    uint32_t method;

    /// @var names - array of names of the individuals in the distance matrix
    /// @details
    /// - if names_src == DMAT_NAMES_SRC_IN_DM_FILE, names is allocated and read from the distance matrix input file
    /// - if names_src == DMAT_NAMES_SRC_IN_VCF_PARS_PTR, names is pointer to the strArray* names in pars which was read from the VCF file, used when metadata is null
    /// - if names_src == DMAT_NAMES_SRC_IN_METADATA_NAMES_PTR, names is pointer to the strArray* names in metadata
    /// - if names_src == DMAT_NAMES_SRC_NONE, names is NULL
    strArray* names;
    uint8_t names_src;

    /// @var matrix - array of n 1D distance matrices of size 'size'
    // matrix[n][size]
    double** matrix;

    /// @var drop[n][size] - 2d array of bool indicators for excluding specific distances from downstream analyses
    /// @note allocated iff args->allow_mispairs==1
    bool** drop;

};

size_t dmat_estimate_max_mem_use(const size_t n, const size_t nInd, const int allow_mispairs, const uint8_t type);
/* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

/* FUNCTION-DECLARATIONS ---------------------------------------------------- */

/// @brief dmat_init - initialize a distance matrix data structure
/// @param n         - number f distance matrices
/// @param nInd      - number of individuals
/// @param type      - type of the distance matrix
/// @param method    - method used to calculate the distances in the matrix
/// @param transform - transformation applied to the distances in the matrix
/// @param names     - array of names of the items in the distance matrix
/// @param names_src - source of the names array
/// @return dmat_t*  - pointer to the newly allocated and initialized distance matrix
dmat_t* dmat_init(const size_t n, const size_t nInd, const uint8_t type, const uint32_t method, uint8_t transform, strArray* names, const uint8_t names_src);

void dmat_destroy(dmat_t* d);

/// @brief dmat_read - read a distance matrix from a distance matrix file
/// @param in_dm_fn           - input distance matrix file name
/// @param required_transform - required transformation
/// @param metadata           - pointer to the metadata data structure
/// @return dmat_t*           - pointer to the distance matrix data structure
dmat_t* dmat_read(const char* in_dm_fn, uint8_t required_transform, metadata_t* metadata);

void dmat_print(dmat_t* dmat, outfile_t* outfile);
void dmat_print_verbose(dmat_t* dmat, outfile_t* outfile);

/// @brief dmat_calculate_distances - calculate the distances using jgtmat and store them in dmat
/// @param jgtmat - pointer to the joint genotype matrix data structure to be used for calculating the distances
/// @param dmat   - pointer to the distance matrix data structure to store the calculated distances
/// @param required_transform - required transformation to be applied to the distances
void dmat_calculate_distances(jgtmat_t* jgtmat, dmat_t* dmat, uint8_t required_transform);

/// @brief new_dmat_pruned - create a new distance matrix data structure with the most dropped pairs removed
/// @details if(dmat->drop[mi]!=NULL); remove the individual with most dropped pairs from the matrix
/// aim: remove as little individuals as possible while ensuring that there is no missing data in the matrix
dmat_t* new_dmat_pruned(dmat_t* dmat);

/// @brief new_dmat_type_to_type - create a new distance matrix data structure with the required type
/// @param dmat          - pointer to the distance matrix data structure
/// @param required_type - required type
/// @return dmat_t*      - pointer to the new distance matrix data structure with the required type
dmat_t* new_dmat_type_to_type(dmat_t* dmat, uint8_t required_type);

/// @brief dmat_type_to_type - change the type of a given distance matrix data structure
/// @param dmat          - pointer to the pointer to the distance matrix data structure
/// @param required_type - required type
/// @return void
void dmat_type_to_type(dmat_t** dmat, uint8_t required_type);

/// @brief dmat_apply_transform - apply a given transformation to the distance matrix in the given distance matrix data structure, if needed
/// @param dmat               - pointer to the distance matrix data structure
/// @param required_transform - required transformation
/// @return void
void dmat_apply_transform(dmat_t* dmat, uint8_t required_transform);

/// @brief new_dmat_matrix_apply_transform - create a new matrix with the required transformation applied 
/// @param dmat               - pointer to the distance matrix data structure
/// @param required_transform - required transformation
/// @return double**          - pointer to the new matrix with the required transformation applied 
double** new_dmat_matrix_apply_transform(dmat_t* dmat, uint8_t required_transform);

/// @brief new_dmat_apply_transform - create a new distance matrix data structure with the required transformation
/// @param dmat               - pointer to the distance matrix data structure
/// @param required_transform - required transformation
/// @return dmat_t*           - pointer to the new distance matrix data structure with the required transformation applied
dmat_t* new_dmat_apply_transform(dmat_t* dmat, uint8_t required_transform);

/// @brief dmat_is_euclidean - check if the distance matrix is euclidean
/// @param dmat - pointer to the distance matrix data structure
/// @param tol  - tolerance for floating point comparison
/// @return bool - true if the distance matrix is euclidean, false otherwise
bool dmat_is_euclidean(dmat_t* dmat, double tol);

/// @brief dmat_has_identity_of_indiscernibles - check if the distance matrix has the identity of indiscernibles property
/// @param dmat - pointer to the distance matrix data structure
/// @return bool - true if the distance matrix has the identity of indiscernibles property, false otherwise
bool dmat_has_identity_of_indiscernibles(dmat_t* dmat);

/* END-OF-FUNCTION-DECLARATIONS --------------------------------------------- */

#endif  // __DMAT_H__