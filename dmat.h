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


/// @brief DMAT_TRANSFORM_* - transformation applied to the distances in the matrix
///
/// NONE: not transformed
#define DMAT_TRANSFORM_NONE 0
/// SQUARE: val^2
#define DMAT_TRANSFORM_SQUARE 1

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
    size_t n;

    /// @var size - size of a each matrix matrix[n]
    /// @details typically nIndCmb
    size_t size;

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

    /// @var has_drop - flag indicating if any distances have been dropped
    /// @details set during distance calculation
    /// @note useful for removing data from downstream analyses
    bool has_drop;

};
/* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

/* FUNCTION-DECLARATIONS ---------------------------------------------------- */
dmat_t* dmat_init(const size_t nInd, const uint8_t type, const uint32_t method, const uint32_t transform, strArray* names, const uint8_t names_src);
void dmat_destroy(dmat_t* d);
dmat_t* dmat_read(const char* in_dm_fn, const uint32_t required_transform, metadata_t* metadata);

void dmat_print(dmat_t* dmat, outfile_t* outfile);
void dmat_print_verbose(dmat_t* dmat, outfile_t* outfile);
void dmat_calculate_distances(jgtmat_t* jgtmat, dmat_t* dmat);


dmat_t* dmat_prune_remove_dropped_distances(dmat_t* dmat);

/* END-OF-FUNCTION-DECLARATIONS --------------------------------------------- */

#endif  // __DMAT_H__