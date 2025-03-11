/**
 * @file    jgtmat.h
 * @brief   header file for jgtmat.cpp
 * @details contains joint genotype matrix data structure (jgtmat_t) and associated functions
 */
#ifndef __JGTMAT_H__
#define __JGTMAT_H__


 /* INCLUDES ----------------------------------------------------------------- */
#include "dataStructs.h"
/* -------------------------------------------------------------------------- */

/* FORWARD DECLARATIONS ----------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* MACROS ------------------------------------------------------------------- */

// JGTMAT9
// calculations based on 3x3 matrix
// 		of pairwise genotype categories
// -------------------------------------
//
// 		MM	Mm	mm
// MM	A	D	G
// Mm	B	E	H
// mm	C	F	I
//
//
// linearized version:
// 		A D G B E H C F I
//      0 1 2 3 4 5 6 7 8
//
// 		MM	Mm	mm
// MM	A=0	D=1	G=2
// Mm	B=3	E=4	H=5
// mm	C=6	F=7	I=8
//

#define JGTMAT_LINEAR_IDXOF_A 0
#define JGTMAT_LINEAR_IDXOF_B 3
#define JGTMAT_LINEAR_IDXOF_C 6
#define JGTMAT_LINEAR_IDXOF_D 1
#define JGTMAT_LINEAR_IDXOF_E 4
#define JGTMAT_LINEAR_IDXOF_F 7
#define JGTMAT_LINEAR_IDXOF_G 2
#define JGTMAT_LINEAR_IDXOF_H 5
#define JGTMAT_LINEAR_IDXOF_I 8

#define JGTMAT_GET_SIJ (A + I + ((B + D + E + F + H) / 2.0))
#define JGTMAT_GET_DIJ (C + G + ((B + D + E + F + H) / 2.0))
#define JGTMAT_GET_FIJ (((2 * C) + (2 * G) - E) / ((2 * C) + (2 * G) + B + D + E + F + H))
#define JGTMAT_GET_IBS0 (C + G)
#define JGTMAT_GET_IBS1 (B + D + F + H)
#define JGTMAT_GET_IBS2 (A + E + I)
#define JGTMAT_GET_R0 ((C + G) / E)
#define JGTMAT_GET_R1 (E / (C + G + B + D + F + H))
#define JGTMAT_GET_KIN ((E - ((2 * C) + (2 * G))) / (B + D + F + H + (2 * E)))

/* END-OF-MACROS ------------------------------------------------------------ */

/* TYPEDEF-STRUCTS ---------------------------------------------------------- */
typedef struct jgtmat_t jgtmat_t;
struct jgtmat_t {

    /// @var n - number of jgtm matrix datasets
    /// @details 1 if no bootstrapping
    ///            args->nBootstraps + 1 (original) if bootstrapping enabled
    size_t n;

    /// @var size - size of each jgtm matrix dataset
    /// @details  typically nIndCmb = (nInd * (nInd - 1)) / 2
    /// @note each dataset contains nIndCmb matrices of size defined by type
    ///       each individual combination in datasets are arranged as LTED
    size_t size;

    /// @var n_gc - number of genotype categories
    /// @details 
    /// - 9 (3x3) for MM Mm mm
    /// - 100 (10x10) for all genotype combinations
    uint8_t n_gc;

    // probabilities
    // pm[n][size][n_gc]
    double*** pm;

    // counts
    // m[n][size][n_gc]
    uint64_t*** m;

    /// @var snsites[n][size] - array of number of sites shared by the two individuals in each matrix
    uint64_t** snsites;

    /// @var drop[n][size] - array of bool indicators for excluding specific matrices from downstream analyses
    /// @note allocated iff args->allow_mispairs==1
    /// @details if gl data, filled in em.cpp. if gt data, filled in vcfReader.cpp
    bool** drop;
};
/* END-OF-TYPEDEF-STRUCTS --------------------------------------------------- */

/* FUNCTION-DECLARATIONS ----------------------------------------------------- */
jgtmat_t* jgtmat_init(const size_t n, const size_t size, const uint8_t n_gc);
void jgtmat_destroy(jgtmat_t* jgtmat);
void jgtmat_print(jgtmat_t* jgtmat, outfile_t* outfile);
size_t jgtmat_estimate_max_mem_use(const size_t n, const size_t size, const uint8_t n_gc);

/// @brief jgtmat_read - read a joint genotype matrix from a joint genotype matrix file
/// @return jgtmat_t* - pointer to the joint genotype matrix data structure
jgtmat_t* jgtmat_read(void);
/* END-OF-FUNCTION-DECLARATIONS ---------------------------------------------- */

#endif // __JGTMAT_H__